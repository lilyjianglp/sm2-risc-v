# SM4 Secure Message Demo

## Purpose

This module adds a small secure communication demo on top of the existing SM2 key exchange implementation. It shows the full application-level flow:

1. Run SM2 key exchange.
2. Verify both sides derive the same shared key.
3. Derive a 16-byte SM4 session key.
4. Encrypt a user-provided message with SM4-CBC.
5. Decrypt the ciphertext and compare it with the original plaintext.

The module is intentionally placed above the arithmetic and RISC-V assembly layers. It does not modify Montgomery multiplication, scalar multiplication, Pornin inversion, or any other low-level optimization code.

## SM2 Integration

The demo uses the top-level SM2 KEX APIs already provided by the project:

- `sm2_kex_initiator_gen_RA`
- `sm2_kex_initiator_compute_key`
- `sm2_kex_responder_compute_key`
- `sm2_kex_initiator_gen_RA_mont`
- `sm2_kex_initiator_compute_key_mont`
- `sm2_kex_responder_compute_key_mont`

`src/sm4_channel.cpp` wraps these calls with:

```cpp
Sm2KexResult run_sm2_key_exchange_demo(bool use_optimized_sm2);
```

The wrapper returns `shared_key_a` and `shared_key_b` so that the SM4 demo can verify `KA == KB` before encryption.

## Key Derivation

The SM4 session key is derived from the SM2 shared key as:

```text
SM4_key = SM3(shared_key || "RV-GM-Secure-SM4-Demo")[0:16]
```

The KDF logic is kept in one function:

```cpp
std::vector<uint8_t> derive_sm4_key(const std::vector<uint8_t>& shared_key);
```

## SM4-CBC Flow

The encryption path is:

1. Generate a 16-byte IV.
2. Apply PKCS#7 padding to the plaintext.
3. Encrypt with SM4-CBC.
4. Print IV and ciphertext in hexadecimal.

The decryption path is:

1. Decrypt with the same SM4 key and IV.
2. Remove PKCS#7 padding.
3. Compare decrypted text with the original input.
4. Print PASS or FAIL.

## Build And Run

From the project `src` directory:

```sh
make secure_message_demo
./secure_message_demo
```

Or run directly:

```sh
make run-secure-message-demo
```

The new target is independent from the existing benchmark targets and does not change `make all`, `make bench-all`, or the low-level assembly source files. It links the existing Montgomery + fp ASM path together with `src/asm/pornin_full_inv.S` for optimized Montgomery-domain inversion.

Menu option 3 invokes:

```sh
make run-pornin-inv-asm USE_FULL_INV_ASM=1
```

Menu option 4 invokes:

```sh
make perf-pornin-asm USE_FULL_INV_ASM=1
```

Run `./secure_message_demo` from the `src` directory so these Makefile commands resolve correctly.

## Example

```text
Please input message:
hello risc-v sm2

[SM2 KEX] Running optimized SM2 key exchange...
[SM2 KEX] Shared key generated successfully.
[SM2 KEX] KA == KB: true

[KDF] Deriving SM4 session key from SM2 shared key...
[SM4] Encrypting message...
[SM4] Ciphertext: 9f8231a0...

[SM4] Decrypting message...
[SM4] Plaintext: hello risc-v sm2

[PASS] Secure message demo success.
```

## Layering Note

This feature is an application-level demonstration module. It only calls the existing upper-level SM2 KEX interface and then performs SM4 message protection. For the demo target, the linked SM2 backend is Montgomery + fp ASM + the existing Pornin full inversion ASM (`src/asm/pornin_full_inv.S`). The module does not edit those low-level files, reimplement SM2, or change the existing baseline and optimized correctness or benchmark paths.
