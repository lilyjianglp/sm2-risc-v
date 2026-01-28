#include "sm2_curve.h"
#include "sm2_scalar.h"
#include <stdio.h>
#include <string.h>
#include <stdint.h>

static void set_k_u32_be(uint8_t k[32], uint32_t x)
{
    memset(k, 0, 32);
    k[31] = (uint8_t)(x);
    k[30] = (uint8_t)(x >> 8);
    k[29] = (uint8_t)(x >> 16);
    k[28] = (uint8_t)(x >> 24);
}

static void fp_to_hex(char out_hex[65], const fp_t* a)
{
    uint8_t buf[32];
    fp_to_bytes(buf, a);
    for (int i = 0; i < 32; i++) {
        static const char* hexd = "0123456789abcdef";
        out_hex[2*i+0] = hexd[(buf[i] >> 4) & 0xF];
        out_hex[2*i+1] = hexd[(buf[i] >> 0) & 0xF];
    }
    out_hex[64] = '\0';
}

static void print_point(const char* tag, const sm2_affine_t* p)
{
    char hx[65], hy[65];
    fp_to_hex(hx, &p->x);
    fp_to_hex(hy, &p->y);
    printf("%s\n", tag);
    printf("  x = %s\n", hx);
    printf("  y = %s\n", hy);
}

static int affine_equal(const sm2_affine_t* a, const sm2_affine_t* b)
{
    uint8_t ax[32], ay[32], bx[32], by[32];
    fp_to_bytes(ax, &a->x);
    fp_to_bytes(ay, &a->y);
    fp_to_bytes(bx, &b->x);
    fp_to_bytes(by, &b->y);
    return (memcmp(ax, bx, 32) == 0) && (memcmp(ay, by, 32) == 0) &&
           (a->infinity == b->infinity);
}

int main(void)
{
    sm2_curve_init_once();

    sm2_affine_t G;
    sm2_get_base_affine(&G);

    /* 你可以多测几个值：1,2,3,4,5,7,8,15,16,31,32... */
    uint32_t tests[] = {1,2,3,4,5,7,8,15,16,31,32212};
    int all_ok = 1;

    for (size_t ti = 0; ti < sizeof(tests)/sizeof(tests[0]); ti++) {
        uint32_t kv = tests[ti];
        uint8_t k[32];
        set_k_u32_be(k, kv);

        sm2_jacobian_t Jnaf, Jladder, Jwin;
        sm2_scalar_mul_naf(&Jnaf, &G, k);
        sm2_scalar_mul_ladder_ct(&Jladder, &G, k);
        sm2_scalar_mul_window_ct_schemeB(&Jwin, &G, k);

        sm2_jacobian_t arr[3] = { Jnaf, Jladder, Jwin };
        sm2_affine_t out[3];
        fp_t prefix_buf[3 + 1];
        sm2_batch_normalize(out, arr, 3, prefix_buf);


        int ok01 = affine_equal(&out[0], &out[1]);
        int ok02 = affine_equal(&out[0], &out[2]);

        printf("k = %u\n", kv);
        if (ok01 && ok02) {
            printf("[PASS] all match\n\n");
        } else {
            printf("[FAIL] scalar multiplication mismatch\n");
            print_point("NAF:", &out[0]);
            print_point("Ladder:", &out[1]);
            print_point("Window w=5 CT:", &out[2]);
            printf("\n");
            all_ok = 0;
        }
    }

    return all_ok ? 0 : 1;
}
