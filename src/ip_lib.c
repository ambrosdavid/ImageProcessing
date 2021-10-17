/*Gruppo: 40, David Ambros, 881443 */
#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"
#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

#define tripleFor(mat,i,j,k) for(i=0;i<mat->h; i++) for( j=0; j<mat->w;j++) for( k=0;k<mat->k;k++)

void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0;i<t->h;i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", get_val(t,i,j,l));
            }
            printf("\n");
        }
        printf("\n");
    }
}


void ip_mat_show_stats(ip_mat * t){
    unsigned int k;

    compute_stats(t);

    for(k=0;k<t->k;k++){
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}

ip_mat * bitmap_to_ip_mat(Bitmap * img){
    unsigned int i=0,j=0;

    unsigned char R,G,B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }

    compute_stats(out);

    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t){

    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            bm_set_pixel(b, j,i, (unsigned char) get_val(t,i,j,0),
                         (unsigned char) get_val(t,i,j,1),
                         (unsigned char) get_val(t,i,j,2));
        }
    }
    return b;
}


float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(float media, float std){

    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float num = cos(2*PI*y2)*sqrt(-2.*log(y1));

    return media + num*std;
}



/**** PARTE 1: TIPO DI DATI ip_mat E MEMORIA ****/

/* Inizializza una ip_mat con dimensioni h w e k. Ogni elemento è inizializzato a v.
 * Inoltre crea un vettore di stats per contenere le statische sui singoli canali.
 * */
void ex(char* ex){
    printf("\n Errore %s.....", ex);
    exit(1);
}

ip_mat * ip_mat_create(unsigned int h, unsigned int w,unsigned  int k, float v){
    unsigned int i,j,kk;

    ip_mat* ipmat = malloc(sizeof(ip_mat));
    if(!ipmat)ex("errore malloc ip_mat_create");

    ipmat->w=w;
    ipmat->h=h;
    ipmat->k=k;
    ipmat->stat = malloc(sizeof(stats) * k);
    if(!ipmat->stat)ex("errore malloc ip_mat_create");

    ipmat->data = malloc(sizeof(float**) * h);
    if(!ipmat->data)ex("errore malloc ip_mat_create");
    for( i=0;i<h;i++){
        ipmat->data[i] = malloc(sizeof(float*) * w);
        if(!ipmat->data[i])ex("errore malloc ip_mat_create");
        for( j=0;j<w;j++){
            ipmat->data[i][j] = malloc(sizeof(float) * k);
            if(!ipmat->data[i][j])ex("errore malloc ip_mat_create");
            for( kk=0;kk<k;kk++)
                set_val(ipmat,i,j,kk,v);
        }
    }

    return ipmat;
}

/* Libera la memoria (data, stat e la struttura) */
void ip_mat_free(ip_mat *a){
    if(a) {
        unsigned int i, j;
        for (i = 0; i < a->h; i++){
            for (j = 0; j < a->w; j++)
                free(a->data[i][j]);
            free(a->data[i]);
        }
        free(a->data);
        free(a->stat);
        free(a);
    }
}


/* Calcola il valore minimo, il massimo e la media per ogni canale
 * e li salva dentro la struttura ip_mat stats
 * */
void compute_stats(ip_mat * t){
    unsigned int i,j,k;

    float massimo[3]={FLT_MAX,FLT_MAX,FLT_MAX};
    float minimo[3]={FLT_MIN,FLT_MIN,FLT_MIN};
    float somma[3]={0,0,0};

    tripleFor(t,i,j,k){
                minimo[k]=min(minimo[k],get_val(t,i,j,k));
                massimo[k]=max(massimo[k],get_val(t,i,j,k));
                somma[k]++;
            }
    for( i=0;i<t->k;i++) {
        t->stat[i].min = minimo[i];
        t->stat[i].max = massimo[i];
        t->stat[i].mean = (somma[i]/(t->h*t->w));
    }
}

/* Inizializza una ip_mat con dimensioni w h e k.
 * Ogni elemento è generato da una gaussiana con media mean e varianza var */
void ip_mat_init_random(ip_mat * t, float mean, float var){
    unsigned int i,j,k;
    tripleFor(t,i,j,k)set_val(t,i,j,k,(get_normal_random(mean,sqrt(var))));
}

/* Crea una copia di una ip_mat e lo restituisce in output */
ip_mat * ip_mat_copy(ip_mat * in){
    unsigned int i,j,k;
    ip_mat* copia=ip_mat_create(in->h,in->w,in->k,0);
    tripleFor(in,i,j,k)set_val(copia,i,j,k,get_val(in,i,j,k));

    for( i=0;i<in->k;i++)(copia->stat)[i]=(in->stat)[i];
    return copia;
}

/* Restituisce una sotto-matrice, ovvero la porzione individuata da:
 * t->data[row_start...row_end][col_start...col_end][0...k]
 * La terza dimensione la riportiamo per intero, stiamo in sostanza prendendo un sottoinsieme
 * delle righe e delle colonne.
 * */
ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end)
{
    unsigned int i,j,k;
    ip_mat * submat;
    if(row_start>(t->h)|| row_end>(t->h) || col_start>(t->w)|| col_end>(t->w)){
        ex("Errore in ip_mat_subset");
    }

    submat=ip_mat_create(row_end-row_start, col_end-col_start,t->k, 0);
    for( i=row_start; i<row_end;i++){
        for( j=col_start;j<col_end;j++){
            for( k=0;k<t->k;k++){
                set_val(submat,i,j,k,get_val(t,i,j,k));
            }
        }
    }
    return submat;
}
/* Concatena due ip_mat su una certa dimensione.
 * Ad esempio:
 * ip_mat_concat(ip_mat * a, ip_mat * b, 0);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h + b.h
 *      out.w = a.w = b.w
 *      out.k = a.k = b.k
 *
 * ip_mat_concat(ip_mat * a, ip_mat * b, 1);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h = b.h
 *      out.w = a.w + b.w
 *      out.k = a.k = b.k
 *
 * ip_mat_concat(ip_mat * a, ip_mat * b, 2);
 *      produrrà un nuovo ip_mat di dimensioni:
 *      out.h = a.h = b.h
 *      out.w = a.w = b.w
 *      out.k = a.k + b.k
 * */
ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dim){
    unsigned int i,j,k,newH,newW,newK;
    ip_mat * mat;
    newH = dim==2 ? a->h + b->h : a->h;
    newW = dim==1 ? a->w + b->w : a->w;
    newK = dim==3 ? a->k + b->k : a->k;
    mat=ip_mat_create(newH, newW,newK, 0);

    for( i=0; i<newH;i++) {
        for ( j = 0; j < newW; j++) {
            for ( k = 0; k < newK; k++) {
                if(i>=a->h || j>=a->w || k>=a->k){
                    set_val(mat,i,j,k,get_val(b,i-(a->h) * (dim==2),j-(a->w)* (dim==1),k-(a->k)*(dim==3)));
                }else set_val(mat,i,j,k,get_val(a,i,j,k));
            }
        }
    }
    return mat;

}

/**** PARTE 1: OPERAZIONI MATEMATICHE FRA IP_MAT ****/
/* Esegue la somma di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output. */
ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b){
    unsigned int i,j,k;
    ip_mat * mat;
    if(a->h!=b->h ||a->w!=b->w || a->k!=b->k ){
        ex("Errore ip_mat_sum--->dimensioni non uguali");
    }
    mat=ip_mat_create(a->h, a->w,a->k, 0);
    tripleFor(a,i,j,k)set_val(mat,i,j,k,get_val(a,i,j,k)+get_val(b,i,j,k));
    return mat;
}

/* Esegue la sottrazione di due ip_mat (tutte le dimensioni devono essere identiche)
 * e la restituisce in output.
 * */
ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b){
    unsigned int i,j,k;
    ip_mat * mat=ip_mat_create(a->h, a->w,a->k, 0);
    tripleFor(a,i,j,k)set_val(mat,i,j,k,get_val(a,i,j,k)-get_val(b,i,j,k));
    return mat;
}

/* Moltiplica un ip_mat per uno scalare c. Si moltiplica c per tutti gli elementi di "a"
 * e si salva il risultato in un nuovo tensore in output. */
ip_mat * ip_mat_mul_scalar(ip_mat *a, float c){
    unsigned int i,j,k;
    ip_mat * mat=ip_mat_create(a->h, a->w,a->k, 0);
    tripleFor(a,i,j,k)set_val(mat,i,j,k,get_val(a,i,j,k)*c);
    return mat;
}

/* Aggiunge ad un ip_mat uno scalare c e lo restituisce in un nuovo tensore in output. */
ip_mat *  ip_mat_add_scalar(ip_mat *a, float c){
    unsigned int i,j,k;
    ip_mat * mat=ip_mat_create(a->h, a->w,a->k, 0);
    tripleFor(a,i,j,k)set_val(mat,i,j,k,get_val(a,i,j,k)+c);
    return mat;
}

/* Calcola la media di due ip_mat a e b e la restituisce in output.*/
ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b){
    unsigned int i,j,k;
    ip_mat * mat=ip_mat_create(a->h, a->w,a->k, 0);
    tripleFor(a,i,j,k)set_val(mat,i,j,k,(get_val(a,i,j,k)+get_val(b,i,j,k))/2);
    return mat;
}

/**** PARTE 2: SEMPLICI OPERAZIONI SU IMMAGINI ****/
/* Converte un'immagine RGB ad una immagine a scala di grigio.
 * Quest'operazione viene fatta calcolando la media per ogni pixel sui 3 canali
 * e creando una nuova immagine avente per valore di un pixel su ogni canale la media appena calcolata.
 * Avremo quindi che tutti i canali saranno uguali.
 * */

ip_mat * ip_mat_to_gray_scale(ip_mat * a){
    unsigned int i,j,k;

    ip_mat * mat=ip_mat_create(a->h, a->w,a->k, 0);
    for( i=0; i<a->h;i++) {
        for ( j = 0; j < a->w; j++) {
            float med=0;
            for( k=0;k<a->k;k++) med +=get_val(a,i,j,k);
            med/=a->k;
            for( k=0;k<a->k;k++)set_val(mat,i,j,k,med);
        }
    }
    return mat;
}

/* Effettua la fusione (combinazione convessa) di due immagini */
ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha){
    unsigned int i,j,k;
    ip_mat * mat;
    alpha=min(alpha,1);
    alpha=max(0,alpha);
    mat=ip_mat_create(a->h, a->w,a->k, 0);
    tripleFor(a,i,j,k)set_val(mat,i,j,k,(alpha*get_val(a,i,j,k))+((1-alpha)*get_val(b,i,j,k)));
    return mat;
}

/* Operazione di brightening: aumenta la luminosità dell'immagine
 * aggiunge ad ogni pixel un certo valore*/
ip_mat * ip_mat_brighten(ip_mat * a, float bright){
    return ip_mat_add_scalar(a,bright);
}

/* Operazione di corruzione con rumore gaussiano:
 * Aggiunge del rumore gaussiano all'immagine, il rumore viene enfatizzato
 * per mezzo della variabile amount.
 * out = a + gauss_noise*amount
 * */
float get_old_normal_random(){
    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    return cos(2*PI*y2)*sqrt(-2.*log(y1));
}
ip_mat * ip_mat_corrupt(ip_mat * a, float amount){
    unsigned int i,j,k;
    ip_mat* ris=ip_mat_copy(a);
    tripleFor(a,i,j,k)set_val(ris,i,j,k,get_val(ris,i,j,k)+(amount*get_old_normal_random()));
    return ris;
}

/**** PARTE 3: CONVOLUZIONE E FILTRI *****/
/* Aggiunge un padding all'immagine. Il padding verticale è pad_h mentre quello
 * orizzontale è pad_w.
 * L'output sarà un'immagine di dimensioni:
 *      out.h = a.h + 2*pad_h;
 *      out.w = a.w + 2*pad_w;
 *      out.k = a.k
 * con valori nulli sui bordi corrispondenti al padding e l'immagine "a" riportata
 * nel centro
 * */
ip_mat * ip_mat_padding(ip_mat * a,unsigned int pad_h,unsigned int pad_w){
    unsigned int i,j,k;
    ip_mat * mat=ip_mat_create(a->h+(pad_h*2), a->w+(pad_w*2),a->k, 0);

    for( i=pad_h; i<a->h+(pad_h);i++) {
        for ( j = pad_w; j <a->w+pad_w; j++) {
            for ( k = 0; k < a->k; k++) {
                set_val(mat,i,j,k,get_val(a,i-pad_h,j-pad_w,k));
            }
        }
    }
    return mat;
}
/* Effettua la convoluzione di un ip_mat "a" con un ip_mat "f".
 * La funzione restituisce un ip_mat delle stesse dimensioni di "a".
 * */
ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f){
    unsigned int i,j,k,l,o;
    ip_mat * mat=ip_mat_padding(a, ((f->h-1)/2), ((f->w-1)/2));
    ip_mat * ris=ip_mat_create(a->h, a->w,a->k,0);
    tripleFor(a,i,j,k){
                float media=0;
                for( l=i; l<i+f->h;l++){
                    for( o=j; o<j+f->w;o++){
                        media+=get_val(f,l-i,o-j,k)*get_val(mat,l,o,k);
                    }
                }
                set_val(ris,i,j,k,media);
            }
    ip_mat_free(mat);
    return ris;
}

ip_mat * to_ip_mat(float *mat,int rows,int cols,int canals) {
    ip_mat *f;
    unsigned int i,j,k;
    f = ip_mat_create(rows, cols, canals, 0.0F);
    tripleFor(f,i,j,k)
                set_val(f,i,j,k,(*((mat+i*cols) + j)));
    return f;
}
/* Crea un filtro di sharpening */
ip_mat * create_sharpen_filter(){
    const float mat[3][3] = {{0,-1,0},{-1,5,-1},{0,-1,0}};
    return to_ip_mat((float *)mat,3,3,3);
}

/* Crea un filtro per rilevare i bordi */
ip_mat * create_edge_filter(){
    const float mat[3][3] = {{-1,-1,-1},{-1,8,-1},{-1,-1,-1}};
    return to_ip_mat((float *)mat,3,3,3);
}

/* Crea un filtro per aggiungere profondità */
ip_mat * create_emboss_filter(){
    const float mat[3][3] = {{-2,-1,0},{-1,1,1},{0,1,2}};
    return to_ip_mat((float *)mat,3,3,3);
}

/* Crea un filtro medio per la rimozione del rumore */
ip_mat * create_average_filter(unsigned int w,unsigned int h,unsigned int k){
    ip_mat * t ;    unsigned int l,j,i;

    t = ip_mat_create(h,w,k,0);
    for(i=0; i<t->h; i++)
        for(j=0; j<t->w; j++)
            for(l=0; l<t->k; l++)
                set_val(t,i,j,l,(1.0/((w*h))));
    return t;
}

/* Crea un filtro gaussiano per la rimozione del rumore */
ip_mat * create_gaussian_filter(unsigned int width,unsigned int height_,unsigned int canals_, float sigma){
    int w=width,h=height_,k=canals_,x,y,i,j;
    double r,s,sommaa,gauss_;
    unsigned int _can;
    ip_mat * output;
    w+=w%2!=0;
    h+=h%2!=0;

    output = ip_mat_create(h+1,w+1,k,0.0F);
    s=2.0*sigma*sigma;
    sommaa=0.0;
    for(x=-(h/2);x<=h/2;x++) {
        for(y=-(w/2);y<=w/2;y++) {
            r=sqrt(x*x+y*y);
            gauss_ = (exp(-(r*r)/s))/(PI*s);
            sommaa+=gauss_;
            for(_can=0;_can<output->k;_can++)set_val(output,x+(h/2),y+(w/2) , _can , gauss_);
        }
    }
    for (i=0;i<h;i++){
        for (j=0;j<w;j++){
            for(_can=0;_can<output->k;_can++)set_val(output,i,j,_can,get_val(output,i,j,_can)/sommaa);
        }
    }
    return output;
}
/*** tirare fuori tutti gli i j nei for**/

ip_mat * spec(ip_mat * a){
    unsigned int i,j,k;
    ip_mat * mat=ip_mat_create(a->h, a->w,a->k, 0);

    for( i=0; i<a->h;i++) {
        for ( j = 0; j < a->w ; j++) {
            for( k=0;k<a->k;k++){
                set_val(mat,i,(a->w)-(1+j),k,get_val(a,i,j,k));
            }
        }
    }
    return mat;
}

ip_mat * mine(ip_mat * a,int amount){
    int jj=0;
    unsigned int i,j,k;
    ip_mat * mat=ip_mat_create(a->h, a->w,a->k, 0);

    for( i=0; i<a->h;i++) {
        for ( j = 0; j < (a->w); j++) {
            if(j%amount==0)jj=j;
            for( k=0;k<(a->k);k++){
                set_val(mat,i,j,k,get_val(a,i,jj,k));
            }
        }
    }
    return mat;
}

/* Effettua una riscalatura dei dati tale che i valori siano in [0,new_max].
 * Utilizzate il metodo compute_stat per ricavarvi il min, max per ogni canale.
 *
 * I valori sono scalati tramite la formula valore-min/(max - min)
 *
 * Si considera ogni indice della terza dimensione indipendente, quindi l'operazione
 * di scalatura va ripetuta per ogni "fetta" della matrice 3D.
 * Successivamente moltiplichiamo per new_max gli elementi della matrice in modo da ottenere un range
 * di valori in [0,new_max].
 * */
void rescale(ip_mat * t, float new_max){
    unsigned int i,j,k;
    compute_stats(t);

    for( k=0;k<t->k;k++)
        for( i=0;i<t->h;i++)
            for( j=0;j<t->w;j++)
                set_val(t,i,j,k,((get_val(t,i,j,k)-((t->stat[k].min)))/((t->stat[k].max)-(t->stat[k].min))*new_max));
}

/* Nell'operazione di clamping i valori <low si convertono in low e i valori >high in high.*/
void clamp(ip_mat * t, float low, float high){
    unsigned int i,j,k;

    tripleFor(t,i,j,k){
                set_val(t,i,j,k,max(low, get_val(t,i,j,k)));
                set_val(t,i,j,k,min(high, get_val(t,i,j,k)));
            }
}


/*
 *
 *
gcc ../bmp.c ../ip_lib.c -c -Wall -lm
gcc ../main_iplib.c ip_lib.o bmp.o -lm

a.exe flower.bmp fullmoon.bmp avg avgg.bmp 1 20 20
a.exe flower.bmp fullmoon.bmp blend blendd.bmp 1 20 0.5
a.exe flower.bmp fullmoon.bmp brighten brightenn.bmp 1 10
a.exe flower.bmp fullmoon.bmp corrupt corruptt.bmp 1 60
a.exe flower.bmp fullmoon.bmp edge edgee.bmp 1
a.exe flower.bmp fullmoon.bmp emboss embosss.bmp 1 2
a.exe flower.bmp fullmoon.bmp gauss gausss.bmp 1 10 2
a.exe flower.bmp fullmoon.bmp gray grayy.bmp 1
a.exe flower.bmp fullmoon.bmp sharp sharpp.bmp 1




 */

