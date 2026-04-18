/* ============================================================
 * sm2_kex_mont.c — Montgomery域KEX实现
 * ============================================================ */
#include "sm2_kex.h"
#include "fp.h"
#include "sm2_scalar.h"
#include <string.h>
#include "fn.h"

/* ============================================================
 *  SM3 + KDF（与sm2_kex.c完全相同的副本）
 * ============================================================ */
typedef struct {
    uint32_t h[8]; uint64_t nbits; uint8_t buf[64]; size_t buf_len;
} sm3_ctx;

static uint32_t rol32(uint32_t x, uint32_t n) { return (x<<n)|(x>>(32-n)); }
static uint32_t load_be32(const uint8_t *p) {
    return ((uint32_t)p[0]<<24)|((uint32_t)p[1]<<16)|((uint32_t)p[2]<<8)|(uint32_t)p[3]; }
static void store_be32(uint8_t *p, uint32_t x) {
    p[0]=(uint8_t)(x>>24);p[1]=(uint8_t)(x>>16);p[2]=(uint8_t)(x>>8);p[3]=(uint8_t)x; }
static uint32_t sm3_p0(uint32_t x){return x^rol32(x,9)^rol32(x,17);}
static uint32_t sm3_p1(uint32_t x){return x^rol32(x,15)^rol32(x,23);}
static uint32_t sm3_ff(uint32_t x,uint32_t y,uint32_t z,int j){
    return(j<16)?(x^y^z):((x&y)|(x&z)|(y&z));}
static uint32_t sm3_gg(uint32_t x,uint32_t y,uint32_t z,int j){
    return(j<16)?(x^y^z):((x&y)|((~x)&z));}
static void sm3_compress(uint32_t H[8],const uint8_t block[64]){
    uint32_t W[68],W1[64];
    for(int i=0;i<16;i++)W[i]=load_be32(block+4*i);
    for(int j=16;j<68;j++){uint32_t x=W[j-16]^W[j-9]^rol32(W[j-3],15);W[j]=sm3_p1(x)^rol32(W[j-13],7)^W[j-6];}
    for(int j=0;j<64;j++)W1[j]=W[j]^W[j+4];
    uint32_t A=H[0],B=H[1],C=H[2],D=H[3],E=H[4],F=H[5],G=H[6],HH=H[7];
    for(int j=0;j<64;j++){
        uint32_t Tj=(j<16)?0x79CC4519u:0x7A879D8Au;
        uint32_t SS1=rol32((rol32(A,12)+E+rol32(Tj,j))&0xFFFFFFFFu,7);
        uint32_t SS2=SS1^rol32(A,12);
        uint32_t TT1=(sm3_ff(A,B,C,j)+D+SS2+W1[j])&0xFFFFFFFFu;
        uint32_t TT2=(sm3_gg(E,F,G,j)+HH+SS1+W[j])&0xFFFFFFFFu;
        D=C;C=rol32(B,9);B=A;A=TT1;HH=G;G=rol32(F,19);F=E;E=sm3_p0(TT2);
    }
    H[0]^=A;H[1]^=B;H[2]^=C;H[3]^=D;H[4]^=E;H[5]^=F;H[6]^=G;H[7]^=HH;
}
static void sm3_init(sm3_ctx *ctx){
    ctx->h[0]=0x7380166F;ctx->h[1]=0x4914B2B9;ctx->h[2]=0x172442D7;ctx->h[3]=0xDA8A0600;
    ctx->h[4]=0xA96F30BC;ctx->h[5]=0x163138AA;ctx->h[6]=0xE38DEE4D;ctx->h[7]=0xB0FB0E4E;
    ctx->nbits=0;ctx->buf_len=0;
}
static void sm3_update(sm3_ctx *ctx,const uint8_t *in,size_t inlen){
    ctx->nbits+=(uint64_t)inlen*8;
    while(inlen>0){
        size_t n=64-ctx->buf_len;if(n>inlen)n=inlen;
        memcpy(ctx->buf+ctx->buf_len,in,n);ctx->buf_len+=n;in+=n;inlen-=n;
        if(ctx->buf_len==64){sm3_compress(ctx->h,ctx->buf);ctx->buf_len=0;}
    }
}
static void sm3_final(sm3_ctx *ctx,uint8_t out[32]){
    uint8_t pad[72];pad[0]=0x80;
    size_t zlen=(ctx->buf_len<56)?(56-ctx->buf_len-1):(64+56-ctx->buf_len-1);
    memset(pad+1,0,zlen);
    uint64_t nbits=ctx->nbits;
    store_be32(pad+1+zlen+0,(uint32_t)(nbits>>32));
    store_be32(pad+1+zlen+4,(uint32_t)(nbits&0xFFFFFFFFu));
    sm3_update(ctx,pad,1+zlen+8);
    for(int i=0;i<8;i++)store_be32(out+4*i,ctx->h[i]);
}
static void sm3_kdf(uint8_t *out,size_t klen,const uint8_t *Z,size_t Zlen){
    uint32_t ct=1;uint8_t buf[32],tmp[4];size_t off=0;
    while(off<klen){
        tmp[0]=(uint8_t)(ct>>24);tmp[1]=(uint8_t)(ct>>16);tmp[2]=(uint8_t)(ct>>8);tmp[3]=(uint8_t)ct;
        sm3_ctx c;sm3_init(&c);sm3_update(&c,Z,Zlen);sm3_update(&c,tmp,4);sm3_final(&c,buf);
        size_t n=(klen-off<32)?(klen-off):32;memcpy(out+off,buf,n);off+=n;ct++;
    }
}

/* ============================================================
 *  曲线参数 + 点合法性（使用 Montgomery 域运算）
 * ============================================================ */
static const uint8_t SM2_A_BYTES[32]={0xFF,0xFF,0xFF,0xFE,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0x00,0x00,0x00,0x00,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFF,0xFC};
static const uint8_t SM2_B_BYTES[32]={0x28,0xE9,0xFA,0x9E,0x9D,0x9F,0x5E,0x34,0x4D,0x5A,0x9E,0x4B,0xCF,0x65,0x09,0xA7,0xF3,0x97,0x89,0xF5,0x15,0xAB,0x8F,0x92,0xDD,0xBC,0xBD,0x41,0x4D,0x94,0x0E,0x93};
static const uint8_t SM2_GX_BYTES[32]={0x32,0xC4,0xAE,0x2C,0x1F,0x19,0x81,0x19,0x5F,0x99,0x04,0x46,0x6A,0x39,0xC9,0x94,0x8F,0xE3,0x0B,0xBF,0xF2,0x66,0x0B,0xE1,0x71,0x5A,0x45,0x89,0x33,0x4C,0x74,0xC7};
static const uint8_t SM2_GY_BYTES[32]={0xBC,0x37,0x36,0xA2,0xF4,0xF6,0x77,0x9C,0x59,0xBD,0xCE,0xE3,0x6B,0x69,0x21,0x53,0xD0,0xA9,0x87,0x7C,0xC6,0x2A,0x47,0x40,0x02,0xDF,0x32,0xE5,0x21,0x39,0xF0,0xA0};

int sm2_point_is_valid(const sm2_affine_t *P) {
    if(!P||sm2_affine_is_infinity(P)) return 0;
    /* 用 Montgomery 域运算验证 y^2 = x^3 + ax + b */
    fp_t a_raw,b_raw; fp_from_bytes(&a_raw,SM2_A_BYTES); fp_from_bytes(&b_raw,SM2_B_BYTES);
    fp_t am,bm,xm,ym; fp_to_mont(&am,&a_raw); fp_to_mont(&bm,&b_raw);
    fp_to_mont(&xm,&P->x); fp_to_mont(&ym,&P->y);
    fp_t y2,x2,x3,ax,rhs;
    fp_mont_sqr(&y2,&ym); fp_mont_sqr(&x2,&xm); fp_mont_mul(&x3,&x2,&xm);
    fp_mont_mul(&ax,&am,&xm); fp_mont_add(&rhs,&x3,&ax); fp_mont_add(&rhs,&rhs,&bm);
    return fp_is_equal(&y2,&rhs); }

static void sm2_u16be(uint8_t out[2],uint16_t x){out[0]=(uint8_t)(x>>8);out[1]=(uint8_t)x;}

int sm2_kex_compute_Z(uint8_t Z[32],const uint8_t *id,size_t id_len,const sm2_affine_t *Ppub){
    if(!Z||!id||!Ppub||!sm2_point_is_valid(Ppub)) return 0;
    uint16_t entl=(uint16_t)(id_len*8);uint8_t entl_be[2];sm2_u16be(entl_be,entl);
    uint8_t a_bytes[32],b_bytes[32];
    memcpy(a_bytes,SM2_A_BYTES,32); memcpy(b_bytes,SM2_B_BYTES,32);
    sm3_ctx c;sm3_init(&c);
    sm3_update(&c,entl_be,2);sm3_update(&c,id,id_len);
    sm3_update(&c,a_bytes,32);sm3_update(&c,b_bytes,32);
    sm3_update(&c,SM2_GX_BYTES,32);sm3_update(&c,SM2_GY_BYTES,32);
    uint8_t x_bytes[32],y_bytes[32];
    fp_to_bytes(x_bytes,&Ppub->x);fp_to_bytes(y_bytes,&Ppub->y);
    sm3_update(&c,x_bytes,32);sm3_update(&c,y_bytes,32);
    sm3_final(&c,Z);return 1;}

static void sm2_kex_x_dash(fn_t *x_dash,const fp_t *x_fp){
    uint8_t xb[32];fp_to_bytes(xb,x_fp);fn_x_dash_w127(x_dash,xb);}
static int sm2_compute_t(uint8_t t_be[32],const uint8_t d_be[32],
    const uint8_t r_be[32],const fn_t *x_dash){
    return fn_compute_t(t_be,d_be,r_be,x_dash);}
static void sm2_kex_compute_inner(uint8_t inner[32],const sm2_affine_t *V,
    const uint8_t ZA[32],const uint8_t ZB[32],
    const sm2_affine_t *RA,const sm2_affine_t *RB){
    uint8_t xV[32],yV[32],x1[32],y1[32],x2[32],y2[32];
    fp_to_bytes(xV,&V->x);fp_to_bytes(yV,&V->y);
    fp_to_bytes(x1,&RA->x);fp_to_bytes(y1,&RA->y);
    fp_to_bytes(x2,&RB->x);fp_to_bytes(y2,&RB->y);
    sm3_ctx c;sm3_init(&c);
    sm3_update(&c,xV,32);sm3_update(&c,ZA,32);sm3_update(&c,ZB,32);
    sm3_update(&c,x1,32);sm3_update(&c,y1,32);sm3_update(&c,x2,32);sm3_update(&c,y2,32);
    sm3_final(&c,inner);}
static void sm2_kex_compute_S(uint8_t Sout[32],uint8_t prefix,
    const sm2_affine_t *V,const uint8_t inner[32]){
    uint8_t yV[32];fp_to_bytes(yV,&V->y);
    sm3_ctx c;sm3_init(&c);
    sm3_update(&c,&prefix,1);sm3_update(&c,yV,32);sm3_update(&c,inner,32);
    sm3_final(&c,Sout);}

/* ============================================================
 *  共享点计算（Montgomery域）
 * ============================================================ */
static int sm2_compute_shared_point_mont_full(
    sm2_affine_t *out_xy, const uint8_t t_self_be[32],
    const fn_t *x_other_dash,
    const sm2_affine_t *R_other_mont, const sm2_affine_t *P_other_mont)
{
    if(!out_xy||!t_self_be||!x_other_dash||!R_other_mont||!P_other_mont) return 0;
    uint8_t xdash_be[32]; fn_to_be(xdash_be,x_other_dash);
    sm2_jacobian_t J1;
    sm2_scalar_mul_window_ct_schemeB_mont_affine_to_jacobian_mont(&J1,R_other_mont,xdash_be);
    sm2_add_ja_mont(&J1,&J1,P_other_mont);
    if(sm2_jacobian_is_infinity(&J1)) return 0;
    sm2_jacobian_t J2;
    sm2_scalar_mul_window_ct_schemeB_mont_jacobian_to_jacobian_mont(&J2,&J1,t_self_be);
    if(sm2_jacobian_is_infinity(&J2)) return 0;
    sm2_jacobian_to_affine_mont(out_xy,&J2);
    return !sm2_affine_is_infinity(out_xy);
}

/* ============================================================
 *  协议：Montgomery域
 * ============================================================ */
int sm2_kex_initiator_gen_RA_mont(sm2_affine_t *RA, const uint8_t rA[32]) {
    if(!RA||!rA) return 0;
    sm2_curve_init_once_mont();
    sm2_affine_t G; sm2_get_base_affine(&G);
    sm2_scalar_mul_window_ct_schemeB_mont_to_affine(RA,&G,rA);
    return sm2_point_is_valid(RA); }

int sm2_kex_initiator_compute_key_mont(uint8_t *K, size_t klen,
    uint8_t S1[32], uint8_t SA[32], const uint8_t *peer_S2,
    const uint8_t dA[32], const sm2_affine_t *PA,
    const uint8_t *idA, size_t idA_len, const uint8_t rA[32],
    const sm2_affine_t *RA, const sm2_affine_t *PB,
    const uint8_t *idB, size_t idB_len, const sm2_affine_t *RB)
{
    if(!K||!S1||!SA||!dA||!PA||!idA||!rA||!RA||!PB||!idB||!RB) return 0;
    if(!sm2_point_is_valid(PA)||!sm2_point_is_valid(PB)||
       !sm2_point_is_valid(RA)||!sm2_point_is_valid(RB)) return 0;
    uint8_t ZA[32],ZB[32];
    if(!sm2_kex_compute_Z(ZA,idA,idA_len,PA)) return 0;
    if(!sm2_kex_compute_Z(ZB,idB,idB_len,PB)) return 0;
    fn_t x1d,x2d; sm2_kex_x_dash(&x1d,&RA->x); sm2_kex_x_dash(&x2d,&RB->x);
    uint8_t tA[32]; if(!sm2_compute_t(tA,dA,rA,&x1d)) return 0;
    sm2_affine_t RB_mont,PB_mont;
    sm2_affine_to_mont(&RB_mont,RB); sm2_affine_to_mont(&PB_mont,PB);
    sm2_affine_t V;
    if(!sm2_compute_shared_point_mont_full(&V,tA,&x2d,&RB_mont,&PB_mont)) return 0;
    uint8_t xV[32],yV[32]; fp_to_bytes(xV,&V.x); fp_to_bytes(yV,&V.y);
    uint8_t buf[128]; memcpy(buf,xV,32); memcpy(buf+32,yV,32);
    memcpy(buf+64,ZA,32); memcpy(buf+96,ZB,32);
    sm3_kdf(K,klen,buf,128);
    uint8_t inner[32]; sm2_kex_compute_inner(inner,&V,ZA,ZB,RA,RB);
    sm2_kex_compute_S(S1,0x02,&V,inner); sm2_kex_compute_S(SA,0x03,&V,inner);
    if(peer_S2&&memcmp(peer_S2,SA,32)!=0) return 0;
    return 1; }

int sm2_kex_responder_compute_key_mont(sm2_affine_t *RB,
    uint8_t *K, size_t klen, uint8_t S2[32], uint8_t SB[32],
    const uint8_t peer_S1[32], const uint8_t dB[32],
    const sm2_affine_t *PB, const uint8_t *idB, size_t idB_len,
    const uint8_t rB[32], const sm2_affine_t *RA,
    const sm2_affine_t *PA, const uint8_t *idA, size_t idA_len)
{
    if(!RB||!K||!S2||!SB||!peer_S1||!dB||!PB||!idB||!rB||!RA||!PA||!idA) return 0;
    if(!sm2_point_is_valid(PA)||!sm2_point_is_valid(PB)||!sm2_point_is_valid(RA)) return 0;
    if(!sm2_kex_initiator_gen_RA_mont(RB,rB)) return 0;
    uint8_t ZA[32],ZB[32];
    if(!sm2_kex_compute_Z(ZA,idA,idA_len,PA)) return 0;
    if(!sm2_kex_compute_Z(ZB,idB,idB_len,PB)) return 0;
    fn_t x1d,x2d; sm2_kex_x_dash(&x1d,&RA->x); sm2_kex_x_dash(&x2d,&RB->x);
    uint8_t tB[32]; if(!sm2_compute_t(tB,dB,rB,&x2d)) return 0;
    sm2_affine_t RA_mont,PA_mont;
    sm2_affine_to_mont(&RA_mont,RA); sm2_affine_to_mont(&PA_mont,PA);
    sm2_affine_t U;
    if(!sm2_compute_shared_point_mont_full(&U,tB,&x1d,&RA_mont,&PA_mont)) return 0;
    uint8_t xU[32],yU[32]; fp_to_bytes(xU,&U.x); fp_to_bytes(yU,&U.y);
    uint8_t buf[128]; memcpy(buf,xU,32); memcpy(buf+32,yU,32);
    memcpy(buf+64,ZA,32); memcpy(buf+96,ZB,32);
    sm3_kdf(K,klen,buf,128);
    uint8_t inner[32]; sm2_kex_compute_inner(inner,&U,ZA,ZB,RA,RB);
    uint8_t expect_S1[32]; sm2_kex_compute_S(expect_S1,0x02,&U,inner);
    if(memcmp(expect_S1,peer_S1,32)!=0) return 0;
    sm2_kex_compute_S(S2,0x03,&U,inner); memcpy(SB,expect_S1,32);
    return 1; }