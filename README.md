> **Important:** This code was developed as part of a Bachelor's Final Degree Project.

***ML-KEM Educational Implementation***

This repository contains an educational implementation of ML-KEM (Module-Lattice-Based Key Encapsulation Mechanism) as specified in FIPS 203.

**Purpose**

The objectives of this project are:

-- To understand the internal structure and design of ML-KEM

-- To study polynomial arithmetic, the Number Theoretic Transform (NTT), compression, and encoding mechanisms

-- To validate the implementation against official Known Answer Tests (KATs)

-- To explore deterministic key generation and encapsulation procedures

This project is intended solely for educational and research purposes.

**Disclaimer**

This implementation is NOT intended for production use.

Specifically, this code:

-- Has not been independently audited

-- Has not been reviewed for side-channel resistance

-- Is not hardened for constant-time execution

-- Does not include production-grade memory protections

-- Has not undergone formal security evaluation

Do NOT Use This Code In

-- Production systems

-- Cryptographic libraries

-- Security-critical environments

-- Real-world deployments

If you require ML-KEM for production use, you should rely on well-reviewed, actively maintained, and professionally audited cryptographic libraries.

**TEST**

This repository includes support for running the official Known Answer Tests (KATs) for ML-KEM (FIPS 203) using the provided .rsp test vector files.

The following code:

-- Reads the values d, z, msg, and ss from each block of the KAT file.

Executes:

* ML_KEMKEYGENINTERN

* ML_KEMENCAPSINTER

* ML_KEMDECAPSINTER

-- Compares the generated shared secret against the expected value (ss).

-- Automatically verifies all blocks in the file.

If all tests pass, the program prints:
ALL <N> BLOCKS OK

**Example**

Place the official KAT file (for example):
kat_MLKEM_768.rsp
in the same directory as the executable.

Compile your project normally.

`Run the executable:
./your_program_name`

The program will process all test blocks sequentially and verify both encapsulation and decapsulation results.

Below is an example main implementation used to run the ML-KEM KAT tests:

```
void hex_to_bytes(const char *hex, uint8_t *out, size_t len)
{
    for (size_t i = 0; i < len; i++) {
        sscanf(hex + 2*i, "%2hhx", &out[i]);
    }
}

int read_kat_block(FILE *f,
                   uint8_t d[32],
                   uint8_t z[32],
                   uint8_t m[32],
                   uint8_t ss[32])
{
    char line[4096];

    int found_d = 0, found_z = 0, found_m = 0, found_ss = 0;

    while (fgets(line, sizeof(line), f)) {

        if (strncmp(line, "d = ", 4) == 0) {
            hex_to_bytes(line + 4, d, 32);
            found_d = 1;
        }
        else if (strncmp(line, "z = ", 4) == 0) {
            hex_to_bytes(line + 4, z, 32);
            found_z = 1;
        }
        else if (strncmp(line, "msg = ", 6) == 0) {
            hex_to_bytes(line + 6, m, 32);
            found_m = 1;
        }
        else if (strncmp(line, "ss = ", 5) == 0) {
            hex_to_bytes(line + 5, ss, 32);
            found_ss = 1;
        }

        if (found_d && found_z && found_m && found_ss)
            return 1;
    }

    return 0;
}

int main() {

    FILE *f = fopen("kat_MLKEM_768.rsp", "r");
    if (!f) {
        printf("Error opening KAT file\n");
        return 1;
    }

    uint8_t d[32], z[32], m[32], ss_ref[32];
    uint8_t ek[384*K+32];
    uint8_t dk[768*K+96];
    uint8_t c[32*(du*K+dv)];
    uint8_t Key[32];
    uint8_t Keydecaps[32];

    int block = 0;

    while (read_kat_block(f, d, z, m, ss_ref)) {

        ML_KEMKEYGENINTERN(ek, dk, d, z);
        ML_KEMENCAPSINTER(Key, c, ek, m);
        ML_KEMDECAPSINTER(Keydecaps, dk, c);

        if (memcmp(Key, ss_ref, 32) != 0) {
            printf("ENCAPS FAIL in block %d\n", block);
            return 1;
        }

        if (memcmp(Keydecaps, ss_ref, 32) != 0) {
            printf("DECAPS FAIL in block %d\n", block);
            return 1;
        }

        printf("Block %d OK\n", block);
        block++;
    }

    fclose(f);

    printf("\nALL %d BLOCKS OK\n", block);
    return 0;
}
```
