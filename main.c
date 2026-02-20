#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h> 
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "randombytes.h"
#include "fips202.h"
#include <x86intrin.h> 
//Los valores d, z y m son declarados de forma constante, si se quiere utilzar las funciones de random quitar los comentarios para d, z y m randombytes
//K,n1,n2 se deben modificar dependiendo para los distintos conjuntos, ML-KEM-512, 768 y 1024
//CPUf depende de la frecuencia del procesador, aunque es solo para medir rendimiento por ciclos
//Parametros criptograficos del esquema ML-KEM
#define N 256
#define Q 3329
#define K 2
#define n1 3
#define n2 2
#define du 10 
#define dv 4
#define di 12
#define CPUf 3700000000UL
//Definicion tipo polinomio: vector de 256 enteros (mod Q)
typedef uint16_t poli[N];
//Valores precalculados para algoritmos 9 y 10
const uint16_t zetas[128] = {
1,1729,2580,3289,2642,630,1897,848,
1062,1919,193,797,2786,3260,569,1746,
296,2447,1339,1476,3046,56,2240,1333,
1426,2094,535,2882,2393,2879,1974,821,
289,331,3253,1756,1197,2304,2277,2055,
650,1977,2513,632,2865,33,1320,1915,
2319,1435,807,452,1438,2868,1534,2402,
2647,2617,1481,648,2474,3110,1227,910,
17,2761,583,2649,1637,723,2288,1100,
1409,2662,3281,233,756,2156,3015,3050,
1703,1651,2789,1789,1847,952,1461,2687,
939,2308,2437,2388,733,2337,268,641,
1584,2298,2037,3220,375,2549,2090,1645,
1063,319,2773,757,2099,561,2466,2594,
2804,1092,403,1026,1143,2150,2775,886,
1722,1212,1874,1029,2110,2935,885,2154
};
//Valores precalculados para algoritmo 11
const int16_t bitrev2[128] = {
17,-17,2761,-2761,583,-583,2649,-2649,
1637,-1637,723,-723,2288,-2288,1100,-1100,
1409,-1409,2662,-2662,3281,-3281,233,-233,
756,-756,2156,-2156,3015,-3015,3050,-3050,
1703,-1703,1651,-1651,2789,-2789,1789,-1789,
1847,-1847,952,-952,1461,-1461,2687,-2687,
939,-939,2308,-2308,2437,-2437,2388,-2388,
733,-733,2337,-2337,268,-268,641,-641,
1584,-1584,2298,-2298,2037,-2037,3220,-3220,
375,-375,2549,-2549,2090,-2090,1645,-1645,
1063,-1063,319,-319,2773,-2773,757,-757,
2099,-2099,561,-561,2466,-2466,2594,-2594,
2804,-2804,1092,-1092,403,-403,1026,-1026,
1143,-1143,2150,-2150,2775,-2775,886,-886,
1722,-1722,1212,-1212,1874,-1874,1029,-1029,
2110,-2110,2935,-2935,885,-885,2154,-2154 
};
//Lectura ciclos
static inline uint64_t read_cycles() {
    return __rdtsc();
}

// =================================================================================
// Suma polinomios
// =================================================================================
void SumPoli(poli out,const uint16_t a[256],const uint16_t b[256]){
    for (int i=0;i<256;i++) {
        out[i]=(a[i]+b[i])%Q;
    }
}
// =================================================================================
// Funcion Traspuesta
// =================================================================================
void Traspuesta(poli A[K][K],poli B[K][K]) {
    for (int i=0;i<K;i++) {
        for (int j=0;j<K;j++) {
            for (int k=0;k<N;k++) {
                A[j][i][k]=B[i][j][k];
            }
        }
    }
}
// =================================================================================
// Funcion Compress
// =================================================================================
void Compress(uint16_t out[],const uint16_t in[],int d) {
    uint16_t mask=(1<<d)-1;
    for (int i=0;i<N;i++) {
        out[i]=(((uint32_t)in[i]<<d)+Q/2)/Q;
        out[i]&=mask;
    }
}

// =================================================================================
// Funcion Decompress
// =================================================================================
void Decompress(uint16_t *a,const uint16_t *b,int d){
    for (int i=0;i<N;i++) {
     a[i]=((uint32_t)b[i]*Q+(1<<(d-1)))>>d; 
    }
    
}
// =================================================================================
// Funcion PRF
// =================================================================================
void PRF(uint8_t out[64*n2], const uint8_t sigma[32],uint8_t NN,int n) {
    uint8_t input[33];
    memcpy(input, sigma, 32);
    input[32]=NN;       
    size_t outlen=64*n;
    shake256(out,outlen,input,33);
}
// =================================================================================
// Algoritmo 3 - BitsToBytes
// =================================================================================
//Convierte array de bits en array de bytes
//256*dl a 32*dl
void BitsToBytes(int dl,uint8_t *B,const uint8_t *bits) {
    size_t BytesL=32*dl;
    memset(B,0,BytesL);
    for (size_t i=0;i<(size_t)256*dl;i++) {
        if (bits[i]) {
            B[i/8]|=(uint8_t)(bits[i]<<(i%8));
        }
    }
}
// =================================================================================
// Algoritmo 4 - BytesToBits
// =================================================================================
//Convierte array de bytes en array de bits
void BytesToBits(int bytelen,uint8_t *b,const uint8_t *B) {
    for (int i=0;i<bytelen;i++) {
        uint8_t C=B[i];
        for (int j=0;j<8;j++) {
            b[i*8+j]=(C>>j)&1u;
        }
    }
}
// =================================================================================
// Algoritmo 5 - ByteEncode
// =================================================================================
//Codifica un arrays de enteros de 𝑑 bits en un array de bytes para 1 ≤ 𝑑 ≤ 12
void ByteEncode(int dl, uint8_t B[32*dl], const uint16_t f[N]) {
    uint8_t bits[256*dl];//Array para almacenar bits(cada uno de los 256 elementos de f se convierte en dl bits)
    for (int i=0;i<256;++i) {
        uint16_t a=f[i];
        for (int j=0;j<dl;++j) {
            bits[i*dl+j]=a&1u;
            a>>=1;
        }
    }
    BitsToBytes(dl,B,bits);
}

// =================================================================================
// Algoritmo 6 - ByteDecode
// =================================================================================
//Decodifica un array de bytes en una array de enteros de 𝑑 bits para 1 ≤ 𝑑 ≤ 12
void ByteDecode(int dl, uint16_t F[N], const uint8_t B[32*dl]) {
    uint8_t bits[256*dl];
    BytesToBits(32*dl, bits, B);
    for (int i = 0; i < N; ++i) {
        uint32_t aux = 0;
        for (int j = 0; j < dl; ++j) {
            aux |= ((uint32_t)bits[i*dl+j])<<j;
        }
        F[i] = aux;
    }
}
// =================================================================================
// Algoritmo 7 - SampleNTT
// =================================================================================
void SampleNTT(uint16_t out[256], const uint8_t input[34]) {
    int j=0;//Indice de coeficientes generados
    size_t max_bytes = 1344;//Suficiente para generar 256 coeficientes (ver Apéndice B)
    uint8_t *buffer = malloc(max_bytes);//Reservar memoria para el buffer
    if (!buffer) {
        fprintf(stderr, "Error: malloc en SampleNTT\n");
        exit(1);
    }
    shake128(buffer,max_bytes,input,34);
    //Extraer coeficientes del buffer hasta completar el polinomio
    size_t offset=0;//Posición actual en el buffer
    while (j<256&&offset+2<max_bytes) {
        //Leer 3 bytes del buffer
        uint8_t C0=buffer[offset];
        uint8_t C1=buffer[offset + 1];
        uint8_t C2=buffer[offset + 2];
        offset+=3;//Avanzar 3 bytes en el buffer
        uint16_t d1=C0+256*(C1%16);//Primer valor (12 bits)
        uint16_t d2=(C1>>4)+16*C2;//Segundo valor (12 bits)
        if(d1<Q){
            out[j]=d1;//Asigna d1 a out[j]
            j=j+1;//Incrementa j después
        }
        if(j<256&&d2<Q) {
            out[j]=d2;
            j=j+1;//Incrementa j después
        }
    }
    free(buffer);//Liberar memoria
}
// =================================================================================
// Funcion generate_matrix
// =================================================================================
void generate_matrix(poli A[K][K], const uint8_t rho[32]) {
    uint8_t input[34];//32 bytes (rho+j+i)
    memcpy(input,rho,32);//Copiar semilla rho
    //Doble bucle para rellenar matriz A[i][j]
    for (int i=0;i<K;i++) {
        for (int j=0;j<K;j++) {
            input[32]=(uint8_t)j;//Byte 33=j(columna)
            input[33]=(uint8_t)i;//Byte 34=i(fila)
            SampleNTT(A[i][j],input);//Generar polinomio
        }
    }
}
// =================================================================================
// Algoritmo 8 - SamplePolyCBD
// =================================================================================
void SamplePolyCBD(uint16_t *r, const uint8_t *B, int n) {
    uint8_t bits[8*64*n];
    BytesToBits(64*n,bits,B);
    for (int i =0;i<256;i++) {
        int a =0,b=0;
        for (int j=0;j<n;j++) {
            a += bits[i*2*n+j];//n bits para a
            b+=bits[i*2*n+n+j];//n bits para b
        }
        int16_t value=a-b;
        r[i] = (value < 0)?value+Q:value;//Asegurar valor pertence Zq
    }
}
// =================================================================================
// Algoritmo 9 - NTT
// =================================================================================
void NTT(uint16_t f[N]) {
    int i=1;
    for (int len=128;len>=2;len/=2) {
        for (int start=0;start<N;start+=2*len) {
            uint16_t zeta=zetas[i];
            i++;
            for (int j=start;j<start+len;j++) {
                //zeta a 32 bits para evitar desbordamiento en el producto
                uint16_t t=(uint16_t)(((uint32_t)zeta*f[j+len])%Q);
                f[j+len]=(f[j]+Q-t)%Q;
                f[j]=(f[j]+t)%Q;
            }
        }
    }
}
// =================================================================================
// Algoritmo 10- NTTi
// =================================================================================
void NTTi(uint16_t f[N]) {
    const uint16_t INV=3303;
    int i=127;
    for (int len=2;len<=128;len*=2) {
        for (int start =0;start<N;start+=2*len) {
            uint16_t zeta = zetas[i--];
            for (int j = start;j<start+len;j++) {
                uint16_t t =f[j];
                f[j]=(t+f[j+len])%Q;
                
                int16_t diff=f[j+len]-t;
                if (diff<0)diff+=Q;

                f[j+len]=(zeta*diff)%Q;
            }
        }
    }
    //Normalizacion
    for (int i=0;i<N;i++) {
        f[i]=(f[i]*INV)%Q;
    }
}
// =================================================================================
// Algoritmo 12- BaseCaseMultiply
// =================================================================================
void BaseCaseMultiply(uint16_t a0,uint16_t a1,uint16_t b0,uint16_t b1,int16_t gammaoriginal,uint16_t *c0,uint16_t *c1){
    //Asegurar que gamma esta en [0,Q-1]
    uint16_t gamma;
    if(gammaoriginal<0) {
        gamma=(uint16_t)(gammaoriginal+Q);
    } else {
        gamma=(uint16_t)gammaoriginal;
    }
    //Usar uint64_t para prevenir desbordamiento
    uint64_t t0=(uint64_t)a0*b0;
    uint64_t t1=(uint64_t)a1*b1*gamma;
    *c0=(uint16_t)((t0+t1)%Q);
    *c1=(uint16_t)(((uint64_t)a0*b1+(uint64_t)a1*b0)%Q);
}
// =================================================================================
// Algoritmo 11- MultiplyNTTs
// =================================================================================
void MultiplyNTTs(uint16_t c[N],const uint16_t a[N],const uint16_t b[N]){
for(int i=0;i<128;i++){
    int16_t bitrev=bitrev2[i];
    BaseCaseMultiply(a[2*i],a[2*i+1],b[2*i],b[2*i+1],bitrev,&c[2*i],&c[2*i+1]);
}
}
// =================================================================================
// Algoritmo 13 - K-PKE.KeyGen
// =================================================================================
void kpke_keygen(uint8_t *ekpke,uint8_t *dkpke,const uint8_t *d) {
    uint8_t buf[33];
    memcpy(buf,d,32);//Copiar semilla d
    buf[32]=K;//Añadir k como byte final
    uint8_t output[64];//Salida de SHA3-512
    sha3_512(output,buf,33);//G(d||k)
    uint8_t rho[32],sigma[32];//Semillas derivadas
    memcpy(rho,output,32);//rho: para generar A
    memcpy(sigma,output+32,32);//sigma: para generar s, e 
    uint8_t NN=0;
    poli A[K][K];//Matriz A pertenece a (Z_q^N)^{K×K}
    generate_matrix(A,rho);//Generar A usando rho
    poli s[K];//Vector secreto s
    int n=n1;
    uint8_t bufPRFs[64*n];
    for (int i=0;i<K;i++) {
    PRF(bufPRFs,sigma,NN,n);
    SamplePolyCBD(s[i],bufPRFs,n);//Obtener s[i]
    NN++;
    }
    poli e[K];
    uint8_t bufPRFe[64*n];
    for (int j=0;j<K;j++) {
    PRF(bufPRFe,sigma,NN,n);
    SamplePolyCBD(e[j],bufPRFe,n);//Obtener e[i]
    NN++;
    }
    //Llamar funcion NTT para shat y ehat
    poli shat[K];
    poli ehat[K];
    for(int i=0;i<K;i++){
    memcpy(shat[i],s[i],N*sizeof(uint16_t));
    NTT(shat[i]);
    memcpy(ehat[i],e[i],N*sizeof(uint16_t));
    NTT(ehat[i]);
    }
    poli that[K];//t en NTT
    poli aux;
    for (int i=0;i<K;i++) {
    memset(that[i],0,N*sizeof(uint16_t));//Inicializar that a 0
    }
    for (int i=0;i<K;i++) {
    memset(aux,0,N*sizeof(uint16_t));//Inicializar suma acumulada para a 0 
    //Multiplicar A[i][j]*s[j] y acumular resultado en aux
    for (int j=0;j<K;j++){
    poli product;
    MultiplyNTTs(product,A[i][j],shat[j]);
    SumPoli(aux,aux,product);
    }
    SumPoli(that[i],aux,ehat[i]);//that[i]=(A*s)[i]+ehat[i]
    }
    //Clave encapsulacion y desencapsulacion
    //Ejecutar K veces ByteEncode y luego introducir rho(semilla A)
    memset(ekpke,0,384*K+32);
    memset(dkpke,0,384*K);
    int dl=di;
    //Generar ekpke
    for(int i=0;i<K;i++){
    ByteEncode(dl,&ekpke[384*i],that[i]);
    }
    //Introducir rho
    for (int i=0;i<32;i++) {
    ekpke[384*K+i]=rho[i];
    }
    //Generar dkpke
    for(int i=0;i<K;i++){
    ByteEncode(dl,&dkpke[384*i],shat[i]);
    }
    //Eliminacion
    memset(rho,0,32);
    memset(sigma,0,32);
    memset(s,0,sizeof(s));
    memset(e,0,sizeof(e));
    memset(shat,0,sizeof(shat));
    memset(ehat,0,sizeof(ehat));
}
// =================================================================================
// Algoritmo 14 - K-PKE.Encrypt
// =================================================================================
void kpke_encrypt(uint8_t c[32*(du*K+dv)],const uint8_t ek[384*K+32],const uint8_t m[32],const uint8_t r[32]){
    printf("\nINICIO PROCESO ENCAPSULACION\n");
    uint8_t NN=0;
    poli that[K];
    int dl=di;
    for(int i=0;i<K;i++){
    ByteDecode(dl,that[i],&ek[384*i]);//Generacion that con ByteDecode
    }
    uint8_t rho[32];
    for(int i=0;i<32;i++) {
    rho[i]=ek[384*K+i];
    }
    poli Ahat[K][K];//Matriz A pertenece a (Z_q^N)^{K×K}
    generate_matrix(Ahat,rho);//Generar A usando rho
    //Generacion y
    poli y[K];//Vector y
    int ny=n1;
    uint8_t bufPRFy[64*ny];//Buffer para salida de PRF
    for (int i=0;i<K;i++) {
    PRF(bufPRFy,r,NN,ny);
    SamplePolyCBD(y[i],bufPRFy,ny);//Obtener y[i]
    NN++;
    }
    //Generacion e1
    poli e1[K];
    int ne=n2;
    uint8_t bufPRFe1[64*ne];//Buffer para salida de PRF
    for (int i=0;i<K;i++) {
    PRF(bufPRFe1,r,NN,ne);//Funcion SHAKE256
    SamplePolyCBD(e1[i],bufPRFe1,ne);//Obtener e1[i]
    NN++;
    }
    //Generacion e2
    uint16_t e2[N];
    uint8_t bufPRFe2[64*ne];//Buffer para salida de PRF
    PRF(bufPRFe2,r,NN,ne);
    SamplePolyCBD(e2,bufPRFe2,ne);//Obtener e2[i]
    //Generacion yhat
    poli yhat[K];
    for(int i=0;i<K;i++){
    memcpy(yhat[i],y[i],N*sizeof(uint16_t));
    NTT(yhat[i]);
    }
    //Generacion u
    poli u[K];
    poli aux;
    //obtener Ahat traspuesta
    poli Ahattraspuesta[K][K];
    Traspuesta(Ahattraspuesta,Ahat);
    //Calcular Ahattraspuesta*yhat+e1
    for (int i=0;i<K;i++) {
    //Inicializar suma acumulada para a 0 
    memset(aux,0,sizeof(poli));
    //Multiplicar Ahattraspuesta[i][j]*yhat[i][j] y acumular resultado en aux
    for (int j=0;j<K;j++){
    poli product;
    MultiplyNTTs(product,Ahattraspuesta[i][j],yhat[j]);
    SumPoli(aux,aux,product);
    }
    NTTi(aux);//Calcular NTT inversa del resultado
    SumPoli(u[i],aux,e1[i]);//Sumar e1[i] al resultado
    }
    //Generacion mu
    int dlmu=1;
    uint16_t s[N];
    uint16_t mu[N];
    ByteDecode(dlmu,s,m);
    Decompress(mu,s,dlmu);
    //Generacion v
    uint16_t v[N];
    memset(v,0,sizeof(v));//Inicializar v a 0
    //Calcular thattraspuesta*yhat+e2+mu
    //Multiplicar thattraspuesta[i]*yhat[i] y acumular resultado en v
    for (int i=0;i<K;i++){
     poli pr;
     MultiplyNTTs(pr,that[i],yhat[i]);
     SumPoli(v,v,pr);
    }
    //Calcular NTT inversa del resultado
    NTTi(v);
    //Sumar e2 y mu al resultado
    SumPoli(v,v,e2);
    SumPoli(v,v,mu);
    //Generacion c1 y c2
    uint16_t x[N];
    uint8_t c1[K*32*du];
    uint8_t c2[32*dv];
    for(int i=0;i<K;i++){
     memset(x,0,sizeof(x));
     Compress(x,u[i],du);
     ByteEncode(du,&c1[i*32*du],x);
    }
    memset(x,0,sizeof(x));
    Compress(x,v,dv);
    ByteEncode(dv,c2,x);
    memcpy(c,c1,K*32*du); 
    memcpy(c+K*32*du,c2,32*dv);
    //Eliminacion
    memset(rho,0,32);
    memset(y,0,sizeof(y));
    memset(yhat,0,sizeof(yhat));
    memset(e1,0,sizeof(e1));
    memset(e2,0,N);
    memset(mu,0,N);

}
// =================================================================================
// Algoritmo 15 - K-PKE.Decrypt
// =================================================================================
void kpke_decrypt(uint8_t m[32],const uint8_t dkpke[384*K],const uint8_t c[32*(du*K+dv)]){
    printf("\nINICIO PROCESO DESENCAPSULACION\n");
    uint8_t c1[32*du*K];
    memcpy(c1,c,32*du*K);//Recuperacion c1 a partir del texto cifrado
    uint8_t c2[32*dv];
    memcpy(c2,&c[32*du*K],32*dv);//Recuperacion c2 a partir del texto cifrado
    int dc=du;
    poli u[K];
    uint16_t x[N];
    for(int i=0;i<K;i++){
        ByteDecode(dc,x,&c1[32*du*i]);
        Decompress(u[i],x,dc);//Recuperacion vector u
    }
    uint16_t v_compressed[N];
    uint16_t v[N];
    ByteDecode(dv,v_compressed,c2);
    Decompress(v,v_compressed,dv);
    poli shat[K];//Vector s en el dominio NTT
    int ds=di;
    for(int i=0;i<K;i++){
       ByteDecode(ds,shat[i],&dkpke[384*i]);//Obtencion del vector s a partir de dkpke
    }
    poli uhat[K];
    for(int i=0;i<K;i++){
        memcpy(uhat[i],u[i],N*sizeof(uint16_t));
        NTT(uhat[i]);
    }
    uint16_t sum[N];
    memset(sum,0,N*sizeof(uint16_t));
    for(int j=0;j<K;j++){
        poli product;
        MultiplyNTTs(product,shat[j],uhat[j]);//shat[j]*uhat[j]
        SumPoli(sum,sum,product);
    }
    //Convertir resultado usando NTTi
    NTTi(sum);
    //Obtener w
    uint16_t w[N]={0};
    for(int i=0;i<N;i++){
     // Se usa int32_t porque el rango de la resta puede ser (-Q, Q)
     // int16_t no garantiza suficiente margen si ocurre un valor negativo
     int32_t temp=(int32_t)v[i]-(int32_t)sum[i];
     if (temp<0)temp+=Q;
     w[i]=(uint16_t)temp;
    }
    //Obtener mensaje
    int dm=1;
    uint16_t aux[N];
    Compress(aux,w,dm);
    ByteEncode(dm,m,aux);
    //Eliminacion
    memset(shat,0,sizeof(shat));
    memset(w,0,N); 
    
}
// =================================================================================
// Algoritmo 16 - ML-KEM.KeyGen_internal
// =================================================================================
void ML_KEMKEYGENINTERN(uint8_t ek[384*K+32],uint8_t dk[768*K+96],const uint8_t d[32],const uint8_t z[32]){
    uint8_t ekpke[384*K+32];//Clave pública interna
    uint8_t dkpke[384*K];//Clave privada interna
    printf("\nINICIO PROCESO DE GENERACION DE CLAVES\n");
    kpke_keygen(ekpke, dkpke, d);
    //Copiar ekpke a ek (son iguales)
    memcpy(ek,ekpke,384*K+32);
    //Calcular hash SHA3-256 de ek
    uint8_t hek[32];
    sha3_256(hek,ek,384*K+32);
    //Emsamblar dk
    memcpy(dk,dkpke,384*K);
    memcpy(dk+384*K,ek,384*K+32);
    memcpy(dk+768*K+32,hek,32);
    memcpy(dk+768*K+64,z,32);
}
// =================================================================================
// Algoritmo 17 - ML-KEM.Encaps_internal
// =================================================================================
void ML_KEMENCAPSINTER(uint8_t Key[32],uint8_t c[32*(du*K+dv)],const uint8_t ek[384*K+32],const uint8_t m[32]){
    //SHA3_256 de ek se guarda en h_ek
    uint8_t h_ek[32];
    sha3_256(h_ek,ek,384*K+32);//H(ek)
    uint8_t KR[64];//32 bytes del mensaje y 32 bytes de H(ek)
    uint8_t mek[64];
    memcpy(KR,m,32);
    memcpy(KR+32,h_ek,32);//m||H(ek)
    sha3_512(mek,KR,64);//G(m||H(ek))
    uint8_t r[32];//valor random r
    memcpy(r,&mek[32],32);//Ultimos 32 bytes a r
    memcpy(Key,mek,32);//Primeros 32 bytes a K
    printf("\nMensaje\n"); 
    for(int i=0;i<32;i++){
    printf("%02X", m[i]);
    } 
    kpke_encrypt(c,ek,m,r);
    printf("\nClave compartida Key generada con ENCAPS\n"); 
    for(int i=0;i<32;i++){
    printf("%02X", Key[i]);
    }
    printf("\nPrimeros 32 bytes del texto cifrado C generado con ENCAPS\n"); 
    for(int i=0;i<32;i++){
    printf("%02X",c[i]);
    }
}
// =================================================================================
// Algoritmo 18 - ML-KEM.Decaps_internal
// =================================================================================
void ML_KEMDECAPSINTER(uint8_t Keydecaps[32],const uint8_t dk[768*K+96],const uint8_t c[32*(du*K+dv)]){
    size_t offset_dkpke=0;//Posicion de dkpke en dk
    size_t offset_ekpke=384*K;//Posicion de ekpke en dk
    size_t offset_h=offset_ekpke+(384*K+32);//Posicion de h en dk
    size_t offset_z=offset_h+32;//Posicion de z en dk
    uint8_t dkpke[384*K];
    memcpy(dkpke,&dk[offset_dkpke],384*K);
    uint8_t ekpke[384*K+32];
    memcpy(ekpke,&dk[offset_ekpke],384*K+32);
    uint8_t h[32];
    memcpy(h,&dk[offset_h],32);
    uint8_t z[32];
    memcpy(z,&dk[offset_z],32);
    uint8_t m[32];
    kpke_decrypt(m,dkpke,c);//Decrypt para descifrar el mensaje recibido y obtener mensaje original
    printf("\nMensaje descifrado\n");   
    for(int i=0;i<32;i++){ 
     printf("%02X", m[i]);
    }
    //G(m|h)
    uint8_t output[64];
    uint8_t input[64];
    memcpy(input,m,32);
    memcpy(&input[32],h,32);
    sha3_512(output,input,64);
    uint8_t Keyde[32];//Secreto compartido generado por Decaps a partir del mensaje recuperado
    uint8_t r[32];
    memcpy(Keyde, output, 32);
    printf("\nClave compartida Key DECAPS G(m|h)\n");   
    for(int i=0;i<32;i++){ 
      printf("%02X", Keyde[i]);
    }
    memcpy(r,&output[32],32);
    //J(z||c)
    uint8_t Key[32];//Secreto compartido generado por rechazo 
    uint8_t zc[32*(du*K+dv)+32];
    memcpy(zc,z,32);
    memcpy(&zc[32],c,32*(du*K+dv));
    size_t outlen=32;
    shake256(Key,outlen,zc,32+32*(du*K+dv));
    printf("\nClave compartida Key DECAPS J(Z||C)\n");   
    for(int i=0;i<32;i++) printf("%02X", Key[i]);
    printf("\n");
    uint8_t chat[32*(du*K+dv)];
    kpke_encrypt(chat,ekpke,m,r);//Llamar encrypt para generar texto cifrado usando mensaje recuperado
    printf("\nPrimeros 32 bytes del texto cifrado chat resultado de funcion ENCRYPT en DECAPS\n"); 
    for(int i=0;i<32;i++) {
    printf("%02X",chat[i]);
    }
    bool a=false;
    for(int i=0;i<32*(du*K+dv);i++){
        if(c[i]!=chat[i]) { 
            a=true; 
            break; 
        }
    }
    if(!a) {
        memcpy(Keydecaps,Keyde,32);
    } else {
        memcpy(Keydecaps,Key,32);
    }
    printf("\nClave compartida Keydecaps SALIDA DECAPS\n");   
    for(int i=0;i<32;i++){
    printf("%02X", Keydecaps[i]);
    }
    printf("\n");
}
// =================================================================================
// Algoritmo 19 - ML-KEM.KeyGen
// =================================================================================
void ML_KEMKEYGEN(uint8_t ek[384*K+32],uint8_t dk[768*K+96]){
    uint8_t d[32] = {
    0xd6, 0x9c, 0xfc, 0x64, 0xf8, 0x4d, 0x4f, 0x33,
    0xe4, 0xc5, 0x4e, 0x16, 0x6b, 0x7f, 0xf9, 0x28,
    0x3a, 0x39, 0x49, 0x86, 0xa5, 0x39, 0xb2, 0x39,
    0x87, 0xa1, 0x0f, 0x39, 0xd2, 0xd9, 0x68, 0x9b
    };
    uint8_t z[32] = {
    0x6d, 0xe6, 0x2e, 0x34, 0x65, 0xa5, 0x5c, 0x9c,
    0x78, 0xa0, 0x7d, 0x26, 0x5b, 0xe8, 0x54, 0x0b,
    0x3e, 0x58, 0xb0, 0x80, 0x1a, 0x12, 0x4d, 0x07,
    0xff, 0x12, 0xb4, 0x38, 0xd5, 0x20, 0x2e, 0xa0
    };
    /* uint8_t d[32];
    uint8_t z[32];
    randombytes(d,32);
    randombytes(z,32);*/
    ML_KEMKEYGENINTERN(ek,dk,d,z);
    printf("\nek(primeros 32 bytes): ");
    for(int i=0;i<32;i++){
      printf("%02X", ek[i]);
    } 
    printf("\ndk(primeros 32 bytes): ");
    for(int i=0; i<32; i++){
      printf("%02X", dk[i]);
    } 
    printf("\nPROCESO GENERACION DE CLAVES FINALIZADO\n");
    printf("\n------------------------------------------------------------------------------\n");
    memset(d,0,32);
    memset(z,0,32);
}
// =================================================================================
// Algoritmo 20 - ML-KEM.Encaps
// =================================================================================
void ML_KEMENCAPS(uint8_t Key[32],uint8_t c[32*(du*K+dv)],uint8_t ek[384*K+32],size_t eksize){
    //Comprobacion ENCAPS
    size_t expected_ek_size=384*K+32;
    if (eksize !=expected_ek_size) { 
        printf("\nERROR: longitud incorrecta de ek. Esperado %zu, recibido %zu\n", 
               expected_ek_size,eksize);
        exit(1);
    }
    uint8_t test[384*K];
    poli that[K];
    for (int i=0;i<K;i++) {
        ByteDecode(di,that[i],&ek[384*i]);
        ByteEncode(di,&test[384*i],that[i]);
    }
    if (memcmp(test,ek,384*K)!=0) {
        printf("\nERROR: ek contiene coeficientes fuera de rango\n");
        exit(1);
    }
    printf("\nComprobacion ENCAPS correcta\n");
    //---------------------------------------------------------------------------------------------------
    uint8_t m[32] = {
    0x01, 0x21, 0xcb, 0x32, 0xac, 0xd1, 0x87, 0x11,
    0x35, 0xcb, 0x34, 0xe2, 0x9c, 0x1a, 0x0e, 0x26,
    0xcc, 0xc0, 0x01, 0xb9, 0x39, 0xea, 0xfa, 0xac,
    0xc2, 0x8f, 0x13, 0xf1, 0x93, 0x8d, 0xbf, 0x91
    };
    /*uint8_t m[32];
    randombytes(m,32);*/
    ML_KEMENCAPSINTER(Key,c,ek,m);
    printf("\nPROCESO DE ENCAPSULACION FINALIZADO\n");
    printf("\n------------------------------------------------------------------------------\n");
}
// =================================================================================
// Algoritmo 21 - ML-KEM.Decaps
// =================================================================================
void ML_KEMDECAPS(uint8_t Keydecaps[32],const uint8_t dk[768*K+96],const uint8_t c[32*(du*K+dv)],size_t dksize,size_t csize){
    //Comprobacion DECAPS 
    size_t expected_dk_size=768*K+96;
    size_t expected_c_size=32*(du*K+dv);
    if (dksize !=expected_dk_size) { 
        printf("\nERROR: longitud incorrecta de dk. Esperado %zu, recibido %zu\n", 
               expected_dk_size,dksize);
        exit(1);
    }
    if (csize !=expected_c_size) { 
        printf("\nERROR: longitud incorrecta de c. Esperado %zu, recibido %zu\n", 
               expected_c_size,csize);
        exit(1);
    }
    uint8_t test[32];
    uint8_t d_hash[384*K+32];
    memcpy(d_hash,&dk[384*K],384*K+32);
    sha3_256(test,d_hash,384*K+32);
    if(memcmp(test,dk+768*K+32,32)!=0) {
        printf("Error: Hash check fallido\n");
        exit(1);
    }
    printf("\nComprobacion DECAPS correcta\n");
    //---------------------------------------------------------------------------------------------------
    ML_KEMDECAPSINTER(Keydecaps,dk,c);
    printf("\nPROCESO DE DESENCAPSULACION FINALIZADO\n");
    printf("\n------------------------------------------------------------------------------\n");
}
// =================================================================================
// Funcion principal
// =================================================================================
int main() {
    uint8_t ek[384*K+32];
    uint8_t dk[768*K+96];
    //KeyGen
    // =================================================================================
    unsigned long long  totalCycles=0;
    uint64_t startKeyGen = read_cycles();
    ML_KEMKEYGEN(ek,dk);
    uint64_t endKeyGen = read_cycles();
    unsigned long long Ke=endKeyGen - startKeyGen;
    double tKe;
    tKe=(double)Ke/CPUf;
    tKe=tKe*1000;
    printf("\nCiclos de CPU usados KeyGen: %llu\n",Ke);
    printf("\nTiempo KeyGen(ms): %f\n",tKe);
    //Encaps
    // =================================================================================
    uint8_t Key[32];
    memset(Key,0,sizeof(Key));
    uint8_t c[32*(du*K+dv)];
    memset(c,0,sizeof(c));
    size_t eksize=sizeof(ek);
    uint64_t startEncaps = read_cycles();
    ML_KEMENCAPS(Key,c,ek,eksize);
    uint64_t endEncaps = read_cycles();
    unsigned long long  E=endEncaps - startEncaps;
    double tE;
    tE=(double)E/CPUf;
    tE=tE*1000;
    printf("\nCiclos de CPU usados Encaps: %llu\n",E);
    printf("\nTiempo Encaps (ms): %f\n",tE);
    //Decaps
    // =================================================================================
    uint8_t Keydecaps[32];
    size_t dksize=sizeof(dk);
    size_t csize=sizeof(c);
    memset(Keydecaps,0,sizeof(Keydecaps));
    uint64_t startDecaps = read_cycles();
    ML_KEMDECAPS(Keydecaps,dk,c,dksize,csize);
    uint64_t endDecaps = read_cycles();
    unsigned long long  D=endDecaps - startDecaps;
    double tD;
    tD=(double)D/CPUf;
    tD=tD*1000;
    printf("\nCiclos de CPU usados Decaps: %llu\n",D);
    printf("\nTiempo Decaps (ms): %f\n",tD);
    totalCycles=Ke+E+D;
    printf("\nCiclos de CPU usados totales: %llu\n",totalCycles);
    return 0;
}