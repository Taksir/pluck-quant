// pcount: count exact occurrences of a fixed set of 25-mers (a probe panel) in reads.
// Only the panel is indexed (hash table + a small bit filter); every read is scanned once.
// Strand: the panel strings are counted as given, on the read strand only (no canonical k-mers),
// which is what the msBWT path and KMC -b count.
//
//   pcount scan  panel.txt threads out.txt  < reads.fastq     stream FASTQ (pipe a decompressor in)
//   pcount pack  store.bin                  < reads.fastq     write a deduplicated 2-bit read store
//   pcount store panel.txt threads out.txt store.bin          scan a read store
//   pcount gz    panel.txt threads out.txt reads.fastq.gz     stream gzip FASTQ, inflated in-process with isa-l (build with -DWITH_ISAL)
#define _GNU_SOURCE
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#ifdef WITH_ISAL
#include <isa-l/igzip_lib.h>
#endif

#define K 25
static const uint64_t KMASK = (1ULL << (2 * K)) - 1;
static int8_t B2[256];

/* ---------- panel index ---------- */
static uint64_t *tkey; static uint32_t *tidx; static int tbits;
static uint64_t *filt; static int fbits;
static uint32_t nuniq, npanel, *line2u;
#define EMPTY (~0ULL)
static inline uint64_t h1(uint64_t x) { return (x * 0x9E3779B97F4A7C15ULL) >> (64 - tbits); }
static inline uint64_t h2(uint64_t x) { return ((x ^ (x >> 29)) * 0xBF58476D1CE4E5B9ULL) >> (64 - fbits); }
static inline int64_t lookup(uint64_t c) {
    uint64_t f = h2(c);
    if (!((filt[f >> 6] >> (f & 63)) & 1)) return -1;
    for (uint64_t m = (1ULL << tbits) - 1, i = h1(c);; i = (i + 1) & m) {
        if (tkey[i] == c) return tidx[i];
        if (tkey[i] == EMPTY) return -1;
    }
}
static int encode(const char *s, uint64_t *c) {
    uint64_t x = 0;
    for (int i = 0; i < K; i++) { int b = B2[(uint8_t)s[i]]; if (b < 0) return 0; x = (x << 2) | b; }
    *c = x; return 1;
}
static void load_panel(const char *fn) {
    FILE *f = fopen(fn, "r"); if (!f) { perror(fn); exit(1); }
    size_t cap = 1 << 20; uint64_t *codes = malloc(cap * 8); char ln[256]; npanel = 0;
    while (fgets(ln, sizeof ln, f)) {
        if (strlen(ln) < K) continue;
        if (npanel == cap) codes = realloc(codes, (cap *= 2) * 8);
        if (!encode(ln, &codes[npanel])) { fprintf(stderr, "bad k-mer %s", ln); exit(1); }
        npanel++;
    }
    fclose(f);
    for (tbits = 4; (1ULL << tbits) < 4ULL * npanel; tbits++);
    for (fbits = 12; (1ULL << fbits) < 32ULL * npanel; fbits++);
    tkey = malloc((1ULL << tbits) * 8); memset(tkey, 0xff, (1ULL << tbits) * 8);
    tidx = malloc((1ULL << tbits) * 4); filt = calloc((1ULL << fbits) / 64 + 1, 8);
    line2u = malloc(npanel * 4); nuniq = 0;
    for (uint32_t j = 0; j < npanel; j++) {
        uint64_t c = codes[j], m = (1ULL << tbits) - 1, i = h1(c);
        while (tkey[i] != EMPTY && tkey[i] != c) i = (i + 1) & m;
        if (tkey[i] == EMPTY) { tkey[i] = c; tidx[i] = nuniq++; uint64_t g = h2(c); filt[g >> 6] |= 1ULL << (g & 63); }
        line2u[j] = tidx[i];
    }
    free(codes);
    fprintf(stderr, "panel %u k-mers (%u unique), table 2^%d, filter 2^%d bits (%.1f MB)\n",
            npanel, nuniq, tbits, fbits, (1ULL << fbits) / 8e6);
}
static inline void scan_seq(const char *s, size_t n, uint32_t *cnt) {
    uint64_t c = 0; int v = 0;
    for (size_t i = 0; i < n; i++) {
        int b = B2[(uint8_t)s[i]];
        if (b < 0) { v = 0; continue; }
        c = ((c << 2) | b) & KMASK;
        if (++v >= K) { int64_t r = lookup(c); if (r >= 0) cnt[r]++; }
    }
}
static void write_counts(const char *fn, uint32_t **per, int nt) {
    for (int t = 1; t < nt; t++) for (uint32_t u = 0; u < nuniq; u++) per[0][u] += per[t][u];
    FILE *o = fopen(fn, "w");
    for (uint32_t j = 0; j < npanel; j++) fprintf(o, "%u\n", per[0][line2u[j]]);
    fclose(o);
}

/* ---------- streaming FASTQ: one reader cuts chunks at record ends, workers scan ---------- */
#define CH (8 << 20)
#define NQ 16
typedef struct { char *b; size_t n; } chunk_t;
static chunk_t q[NQ]; static int qh, qt, qn, done_reading;
static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t cv_put = PTHREAD_COND_INITIALIZER, cv_get = PTHREAD_COND_INITIALIZER;
static uint64_t nreads_total;
static void *worker(void *arg) {
    uint32_t *cnt = arg; uint64_t nr = 0;
    for (;;) {
        pthread_mutex_lock(&mu);
        while (!qn && !done_reading) pthread_cond_wait(&cv_get, &mu);
        if (!qn) { pthread_mutex_unlock(&mu); break; }
        chunk_t c = q[qh]; qh = (qh + 1) % NQ; qn--;
        pthread_cond_signal(&cv_put); pthread_mutex_unlock(&mu);
        char *p = c.b, *e = c.b + c.n; int l = 0;
        while (p < e) {
            char *nl = memchr(p, '\n', e - p); if (!nl) nl = e;
            if ((l & 3) == 1) { scan_seq(p, nl - p, cnt); nr++; }
            l++; p = nl + 1;
        }
        free(c.b);
    }
    __atomic_add_fetch(&nreads_total, nr, __ATOMIC_RELAXED);
    return 0;
}
static void push(char *b, size_t n) {
    pthread_mutex_lock(&mu);
    while (qn == NQ) pthread_cond_wait(&cv_put, &mu);
    q[qt] = (chunk_t){b, n}; qt = (qt + 1) % NQ; qn++;
    pthread_cond_signal(&cv_get); pthread_mutex_unlock(&mu);
}
#ifdef WITH_ISAL
static struct inflate_state ist; static const uint8_t *gzp, *gze;
static ssize_t gz_read(char *dst, size_t n) {   // inflate straight into the chunk buffer (handles multi-member gzip)
    if (!n) return 0;
    for (;;) {
        ist.next_out = (uint8_t *)dst; ist.avail_out = n;
        if (isal_inflate(&ist) != ISAL_DECOMP_OK) { fprintf(stderr, "inflate error\n"); exit(1); }
        ssize_t got = n - ist.avail_out;
        if (ist.block_state == ISAL_BLOCK_FINISH && ist.avail_in > 8) {
            const uint8_t *nx = ist.next_in; uint32_t av = ist.avail_in;
            isal_inflate_reset(&ist); ist.crc_flag = ISAL_GZIP; ist.next_in = (uint8_t *)nx; ist.avail_in = av;
        }
        if (got || ist.avail_in == 0 || ist.block_state == ISAL_BLOCK_FINISH) return got;
    }
}
static int use_gz;
#endif
static void stream(int fd, int nt, const char *out) {
    pthread_t th[64]; uint32_t *per[64];
    for (int t = 0; t < nt; t++) { per[t] = calloc(nuniq, 4); pthread_create(&th[t], 0, worker, per[t]); }
    char *buf = malloc(CH); size_t have = 0; uint64_t line = 0;
    for (;;) {
        ssize_t r;
#ifdef WITH_ISAL
        if (use_gz) r = gz_read(buf + have, CH - have); else
#endif
        r = read(fd, buf + have, CH - have);
        if (r < 0) { perror("read"); exit(1); }
        have += r;
        if (r == 0 || have == CH) {
            // cut after the last complete 4-line record
            size_t cut = 0; uint64_t l = line;
            for (char *p = buf, *e = buf + have; p < e;) {
                char *nl = memchr(p, '\n', e - p); if (!nl) break;
                p = nl + 1; if ((++l & 3) == 0) cut = p - buf;
            }
            if (r == 0 && have > cut) cut = have;  // tail without final newline
            char *nb = malloc(CH); memcpy(nb, buf + cut, have - cut);
            if (cut) push(buf, cut); else free(buf);
            buf = nb; have -= cut;
            if (r == 0) break;
        }
    }
    free(buf);
    pthread_mutex_lock(&mu); done_reading = 1; pthread_cond_broadcast(&cv_get); pthread_mutex_unlock(&mu);
    for (int t = 0; t < nt; t++) pthread_join(th[t], 0);
    fprintf(stderr, "reads %lu\n", (unsigned long)nreads_total);
    write_counts(out, per, nt);
}

/* ---------- read store: unique reads (len <= 64, no N) in 2-bit, plus a count; the rest as text ---------- */
typedef struct { uint64_t hi, lo; uint32_t n; } ent_t;  // hi carries the length in its top byte
static ent_t *ht; static uint64_t hcap, hn;
static inline uint64_t hh(uint64_t hi, uint64_t lo) { return ((hi * 0x9E3779B97F4A7C15ULL) ^ (lo * 0xC2B2AE3D27D4EB4FULL)) * 0x165667B19E3779F9ULL; }
static void hins(uint64_t hi, uint64_t lo, uint32_t n) {
    uint64_t m = hcap - 1, i = hh(hi, lo) & m;
    while (ht[i].n && !(ht[i].hi == hi && ht[i].lo == lo)) i = (i + 1) & m;
    if (!ht[i].n) { ht[i].hi = hi; ht[i].lo = lo; hn++; }
    ht[i].n += n;
}
static void hgrow(void) {
    ent_t *o = ht; uint64_t oc = hcap; hcap *= 2; ht = calloc(hcap, sizeof *ht); hn = 0;
    for (uint64_t i = 0; i < oc; i++) if (o[i].n) hins(o[i].hi, o[i].lo, o[i].n);
    free(o);
}
static void pack(const char *out) {
    hcap = 1 << 24; ht = calloc(hcap, sizeof *ht);
    FILE *in = stdin; static char l1[4096], l2[4096], l3[4096], l4[4096];
    char sidefn[4096]; snprintf(sidefn, sizeof sidefn, "%s.side", out);
    FILE *side = fopen(sidefn, "w"); uint64_t nr = 0, nside = 0; int lens[65] = {0};
    while (fgets(l1, sizeof l1, in) && fgets(l2, sizeof l2, in) && fgets(l3, sizeof l3, in) && fgets(l4, sizeof l4, in)) {
        size_t n = strcspn(l2, "\r\n"); nr++;
        uint64_t hi = 0, lo = 0; int ok = n <= 64;
        for (size_t i = 0; ok && i < n; i++) {
            int b = B2[(uint8_t)l2[i]]; if (b < 0) { ok = 0; break; }
            hi = (hi << 2) | (lo >> 62); lo = (lo << 2) | b;
        }
        if (!ok) { fwrite(l2, 1, n, side); fputc('\n', side); nside++; continue; }
        hi |= (uint64_t)n << 56; lens[n] = 1;  // 2-bit payload of <=64 bases fits 128 bits; length goes in hi's top byte only when n<=60
        if (n > 60) { fwrite(l2, 1, n, side); fputc('\n', side); nside++; continue; }
        if (2 * (hn + 1) > hcap) hgrow();
        hins(hi, lo, 1);
    }
    fclose(side);
    // groups by length: [u32 len][u64 n][n * ceil(len/4) bytes, base i at byte i/4 bits 6-2(i%4)][n * u32 count]
    FILE *o = fopen(out, "wb"); uint32_t magic = 0x50434B31; fwrite(&magic, 4, 1, o);
    for (uint32_t L = 1; L <= 60; L++) {
        if (!lens[L]) continue;
        uint64_t n = 0; for (uint64_t i = 0; i < hcap; i++) if (ht[i].n && (ht[i].hi >> 56) == L) n++;
        uint32_t B = (L + 3) / 4; fwrite(&L, 4, 1, o); fwrite(&n, 8, 1, o);
        uint8_t rec[16];
        for (uint64_t i = 0; i < hcap; i++) if (ht[i].n && (ht[i].hi >> 56) == L) {
            memset(rec, 0, B);
            for (uint32_t j = 0; j < L; j++) {
                uint32_t sh = 2 * (L - 1 - j);
                uint64_t b = sh >= 64 ? (ht[i].hi >> (sh - 64)) & 3 : (ht[i].lo >> sh) & 3;
                rec[j >> 2] |= b << (6 - 2 * (j & 3));
            }
            fwrite(rec, 1, B, o);
        }
        for (uint64_t i = 0; i < hcap; i++) if (ht[i].n && (ht[i].hi >> 56) == L) fwrite(&ht[i].n, 4, 1, o);
    }
    fclose(o);
    fprintf(stderr, "reads %lu, unique packed %lu (%.1f%%), side %lu\n", (unsigned long)nr, (unsigned long)hn,
            100.0 * hn / nr, (unsigned long)nside);
}
typedef struct { const uint8_t *s; const uint32_t *c; uint32_t L, B; uint64_t a, b; uint32_t *cnt; } job_t;
static void *store_worker(void *arg) {
    job_t *j = arg;
    for (uint64_t r = j->a; r < j->b; r++) {
        const uint8_t *p = j->s + r * j->B; uint32_t w = j->c[r]; uint64_t c = 0;
        for (uint32_t i = 0; i < j->L; i++) {
            c = ((c << 2) | ((p[i >> 2] >> (6 - 2 * (i & 3))) & 3)) & KMASK;
            if (i + 1 >= K) { int64_t x = lookup(c); if (x >= 0) j->cnt[x] += w; }
        }
    }
    return 0;
}
static void scan_store(const char *fn, int nt, const char *out) {
    int fd = open(fn, O_RDONLY); struct stat st; fstat(fd, &st);
    const uint8_t *m = mmap(0, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0), *p = m + 4, *e = m + st.st_size;
    uint32_t *per[64]; for (int t = 0; t < nt; t++) per[t] = calloc(nuniq, 4);
    while (p < e) {
        uint32_t L; uint64_t n; memcpy(&L, p, 4); memcpy(&n, p + 4, 8); p += 12;
        uint32_t B = (L + 3) / 4; const uint8_t *s = p; const uint32_t *c = (const uint32_t *)(p + n * B);
        pthread_t th[64]; job_t jb[64];
        for (int t = 0; t < nt; t++) {
            jb[t] = (job_t){s, c, L, B, n * t / nt, n * (t + 1) / nt, per[t]};
            pthread_create(&th[t], 0, store_worker, &jb[t]);
        }
        for (int t = 0; t < nt; t++) pthread_join(th[t], 0);
        p += n * B + n * 4;
    }
    // side reads (N-containing or longer than 60), as plain lines
    char sidefn[4096]; snprintf(sidefn, sizeof sidefn, "%s.side", fn);
    FILE *f = fopen(sidefn, "r"); char ln[4096];
    if (f) { while (fgets(ln, sizeof ln, f)) scan_seq(ln, strcspn(ln, "\n"), per[0]); fclose(f); }
    write_counts(out, per, nt);
}

int main(int argc, char **argv) {
    memset(B2, -1, sizeof B2); B2['A'] = 0; B2['C'] = 1; B2['G'] = 2; B2['T'] = 3;
    if (argc >= 5 && !strcmp(argv[1], "scan")) { load_panel(argv[2]); stream(0, atoi(argv[3]), argv[4]); }
#ifdef WITH_ISAL
    else if (argc >= 6 && !strcmp(argv[1], "gz")) {
        int fd = open(argv[5], O_RDONLY); struct stat st; fstat(fd, &st);
        gzp = mmap(0, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0); gze = gzp + st.st_size;
        isal_inflate_init(&ist); ist.crc_flag = ISAL_GZIP; ist.next_in = (uint8_t *)gzp; ist.avail_in = st.st_size; use_gz = 1;
        load_panel(argv[2]); stream(-1, atoi(argv[3]), argv[4]);
    }
#endif
    else if (argc >= 3 && !strcmp(argv[1], "pack")) pack(argv[2]);
    else if (argc >= 6 && !strcmp(argv[1], "store")) { load_panel(argv[2]); scan_store(argv[5], atoi(argv[3]), argv[4]); }
    else { fprintf(stderr, "usage: see header\n"); return 1; }
    return 0;
}
