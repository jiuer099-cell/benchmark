#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <openssl/sha.h>

/* Hash the complete SAM record stream with only the added RG field removed.
 * Every other field, field order, record order and newline is retained. */
int main(int argc, char **argv) {
    if (argc != 2) { fprintf(stderr, "usage: rg_record_audit expected_rg < records.sam\n"); return 2; }
    SHA256_CTX sha; SHA256_Init(&sha);
    char *line = NULL; size_t cap = 0; ssize_t len;
    uint64_t records = 0;
    while ((len = getline(&line, &cap, stdin)) >= 0) {
        if (!len || line[len-1] != '\n' || line[0] == '@') goto bad;
        size_t start = 0; unsigned field = 0, rg = 0, kept = 0;
        for (size_t end = 0; end < (size_t)len; ++end) {
            if (line[end] != '\t' && line[end] != '\n') continue;
            ++field;
            if (field > 11 && end-start >= 5 && !memcmp(line+start,"RG:Z:",5)) {
                if (++rg != 1 || end-start-5 != strlen(argv[1]) || memcmp(line+start+5,argv[1],end-start-5)) goto bad;
            } else {
                if (kept++) SHA256_Update(&sha,"\t",1);
                SHA256_Update(&sha,line+start,end-start);
            }
            start = end+1;
        }
        if (field < 12 || rg != 1) goto bad;
        SHA256_Update(&sha,"\n",1); ++records;
    }
    if (ferror(stdin) || !records) goto bad;
    unsigned char digest[SHA256_DIGEST_LENGTH]; SHA256_Final(digest,&sha);
    printf("{\"records\":%llu,\"every_read_has_expected_rg\":true,\"normalized_records_sha256\":\"",(unsigned long long)records);
    for (unsigned i=0;i<SHA256_DIGEST_LENGTH;++i) printf("%02x",digest[i]);
    puts("\"}"); free(line); return 0;
bad:
    fprintf(stderr,"Invalid or unexpected RG/SAM stream at record %llu\n",(unsigned long long)records+1);
    free(line); return 1;
}
