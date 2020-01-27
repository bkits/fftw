#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <pthread.h>
#include <time.h>
#define NUM_THREADS 2

struct thread_data{fftw_plan plan;};// ena struct pou periexei ta arguments pou 8elw na perasw sto thread

struct thread_data thread_data_array[NUM_THREADS];// ftiaxnw ena pinaka oste na mporei to ka8e thread na parei diaforetika arguments

/*
void mult( fftw_complex *c,fftw_complex *a, fftw_complex *b,int n0,int n1,int n2 ) 
 {//polaplasiasmos pinakwn complex ari8mwn c[n0][n2]=a[n0][n1]*b[n1][n2] 
    int i, j, k;
    for (i=0; i<n0; ++i)
       {for (j=0; j<n2; ++j) 
        { c[i*n2+j] = 0.0;
          for (k=0; k<n1; ++k)
          {c[i*n2+j] =c[i*n2+j]+a[i*n1+k] * b[k*n2+j];}}} //oi pinakes einai apo malloc opote diabazontai san monodiastatoi
                                                          //opou px to stoixeio c[i][j] tou pinaka c[n0][n2] einai to c[i*n2+j]
   }

*/

void mult( fftw_complex *c,fftw_complex *a, fftw_complex *b,int n0,int n1) 
 {//polaplasiasmos simoon twn pinakwn
    int i;
    for (i=0; i<(n0*n1); i++)
       {
        c[i] =a[i] * b[i];
       } 
   }



void insert(fftw_complex *a,fftw_complex *b,int n0,int n1,int n11,int n2,int lim,int n)// ena insert opou sto pinaka a[n0][n1] kai b[n1][n2] 
{int i,j;time_t t;srand((unsigned)time(&t)); //ton gemizei me tuxea migadika stoixeia 
for(i=0;i<n0;i++)
 {for(j=0;j<n1;j++)
   {a[i*n+j]=rand()%lim+rand()%lim*I;
   }
 }
for(i=0;i<n11;i++)
 {for(j=0;j<n2;j++)
   {b[i*n+j]=rand()%lim+rand()%lim*I;
   }
 }

}

void insert_int(fftw_complex *a,fftw_complex *b,int n0,int n1,int n11,int n2,int lim,int n)// ena insert opou sto pinaka a[n0][n1] kai b[n1][n2] 
{int i,j;time_t t;srand((unsigned)time(&t)); //ton gemizei me tuxea akeraia stoixeia 
for(i=0;i<n0;i++)
 {for(j=0;j<n1;j++) 
   {a[i*n+j]=rand()%lim+0*I;
   }
 }
for(i=0;i<n11;i++)
 {for(j=0;j<n2;j++)
   {b[i*n+j]=rand()%lim+0*I;
   }
 }
printf("\n");
}

void insert_real(fftw_complex *a,fftw_complex *b,int n0,int n1,int n11,int n2,int lim,int n)// ena insert opou sto pinaka a[n0][n1] kai b[n1][n2] 
{int i,j;time_t t;srand((unsigned)time(&t)); //ton gemizei me tuxea pragmatika stoixeia 
lim--;
for(i=0;i<n0;i++)
 {for(j=0;j<n1;j++) 
   {a[i*n+j]=(rand()%lim+(rand()%999)/1000.0)+0*I;
   }
 }
for(i=0;i<n11;i++)
 {for(j=0;j<n2;j++)
   {b[i*n+j]=(rand()%lim+(rand()%999)/1000.0)+0*I;
   }
 }
printf("\n");
}


void print(fftw_complex *a,int n0,int n1,int n,FILE *f)// ektuposi pinaka a[n0][n1] me morfi (pragmatiko meros,fantastiko meros)
{int i,j;
 
for(i=0;i<n0;i++)
 {for(j=0;j<n1;j++)
   {fprintf(f,"%.3lf,%.3lf ",creal(a[i*n+j]),cimag(a[i*n+j]));
   }fprintf(f,"\n");
 }
}

void print_real(fftw_complex *a,int n0,int n1,int n,FILE *f)// ektuposi pinaka a[n0][n1]
{int i,j;
 
for(i=0;i<n0;i++)
 {for(j=0;j<n1;j++)
   {fprintf(f,"%.3lf ",creal(a[i*n+j]));
   }fprintf(f,"\n");
 }
}

void normalize_2D(fftw_complex *y,int n0,int n2)// mia sunartisi pou sto antistrofo metasximatismo diairei ola ta stoixeia me 
{int i;for(i=0;i<(n0*n2);i++){y[i]=y[i]/(n0*n2);}} //ton ari8mo ton bimatwn tou metasximatismou (n0*n2 epeidi eimaste se
                                                  //disdiastato) mias kai den to kanei apo moni tis i sunartisi tis biblio8ikis  

int plithos(FILE *f)// mia sunartisi pou metraei ola ta stoixeia tou pragmatikou pinaka
{int pli8os;double b;
while(!feof(f))
{fscanf(f,"%lf ",&b);
pli8os++;
}return pli8os;
}

int plithos_complex(FILE *f)// mia sunartisi pou metraei ola ta stoixeia tou complex pinaka
{int pli8os;double im,re;
while(!feof(f))
{fscanf(f,"%lf,%lf ",&re,&im);
pli8os++;
}return pli8os;
}

int error_real(FILE *f)//elenxos dedomenwn
{
double b;
while(!feof(f))
{if (fscanf(f,"%lf ",&b)!=1)
{return 1;}
}
return 0;
}

int error_complex(FILE *f)// elenxos dedomenwn gia ta arxeia complex
{double im,re;
while(!feof(f))
{if(fscanf(f,"%lf,%lf ",&re,&im)!=2)
{return 1;}

}return 0;
}




void insert_file(FILE *f,fftw_complex *a,int n0,int n1,int n)
{int b,c;double d;
for(b=0;b<n0;b++)
{
 for(c=0;c<n1;c++)
 {fscanf(f,"%lf ",&d);a[b*n+c]=d;}
}
}

void insert_file_complex(FILE *f,fftw_complex *a,int n0,int n1,int n)
{int b,c;double re,im;
for(b=0;b<n0;b++)
{
 for(c=0;c<n1;c++)
{fscanf(f,"%lf,%lf ",&re,&im);a[b*n+c]=re+im*I;}
}
}

void zero(fftw_complex *a, int n0,int n1,int m,int n) //gemisma me midenika ta upoloipa stoixeia tou pinaka
{int i,j;
 for(i=0;i<n0;i++)
 {
  for (j=n1;j<n;j++)
   {a[i*n+j]=0+0*I;}
 }
 for (i=(n0*n);i<(m*n);i++) {a[i]=0+0*I;}
}

int input_dim() //diabasma kai epali8eusi enos akeraiou apo to xristi
{int a1,a2,r;
a2=scanf("%d",&a1);getchar();
while(a2!=1) 
{printf("prepei na einai akeraios ari8mos\n");
a2=scanf(" %d",&a1);while((r = getchar()) != '\n');
}
while (a1<=0)
{printf("prepei na einai 8etikos ari8mos\n");
 a1=input_dim();
}

return a1;
}



void *execute(void *threadarg){//i sunartisi pou ekteloun ta threads
struct thread_data *data;
data=(struct thread_data *) threadarg;
fftw_execute(data->plan);// i sunartisi tis biblio8ikis fftw3.h pou ektelei to plano (plano 8a eksgi8ei parakato)
pthread_exit(NULL);}// eksodos tou thread to bazoume gia asfaleia

int main(int argc,char *argv[])
{time_t start1,end;FILE *f;
pthread_t threads[NUM_THREADS];//orizoume ta threads
pthread_attr_t attr;//orizoume oti ta threads 8a exoun kapia attributes
int n0,n1,n2,n3,rc,lim,error,n1_aux,lim1,m,n;
char c,choise='*';
fftw_complex *a,*A,*b,*B,*y,*Y;//edw einai oi pinakes pou 8a ftiaksoume (pinakas--F-->PINAKAS parakatw 8a dilonei metasximatismo)
fftw_plan plan_a,plan_b,plan_y;//edw orizoume ta plana
void *status;//ena orisma pou xreiazete i sunartisi pthread_join tis pthead.h
pthread_attr_init(&attr);
pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);//ftiaxnoume to attribute twn thread na einai joinable diladi to 
                                                            //programa na proxorisei afou teliwsoun prota ta threads

while(choise!='q')
{printf("(a) dimiourgia input files\n");
 printf("(r) random complex\n");
 printf("(s) random integer\n");
 printf("(f) random real\n");
 printf("(p) diabasma apo arxeia real\n");
 printf("(c) diabasma apo arxeia complex\n");
printf("(q) quit\n");
 scanf(" %c",&choise);while((lim1 = getchar()) != '\n');
// getchar();
if(choise=='a')
{printf("\n");
if((f=fopen("input_a","w"))==NULL){printf("ERROR OPENING FILE");exit(1);}
fclose(f);printf("dimiourgia input_a\n");
if((f=fopen("input_a_complex","w"))==NULL){printf("ERROR OPENING FILE");exit(1);}
fclose(f);printf("dimiourgia input_a_complex\n");
if((f=fopen("input_b","w"))==NULL){printf("ERROR OPENING FILE");exit(1);}
fclose(f);printf("dimiourgia input_b\n");
if((f=fopen("input_b_complex","w"))==NULL){printf("ERROR OPENING FILE");exit(1);}
fclose(f);printf("dimiourgia input_b_complex\n\n");
}
if((choise=='r')||(choise=='s')||(choise=='p')||(choise=='c')||(choise=='f'))
{
if (choise=='r')
// EISODOS RAND complex
{
printf("\n dose diastaseis n0xn1 n2xn3\n");
printf(" n0=");n0=input_dim();
printf(" n1=");n1=input_dim();
printf(" n2=");n2=input_dim();
printf(" n3=");n3=input_dim();
printf(" dose orio tuxaion ari8mon ");
lim=input_dim();

m=n0+n2-1;
n=n1+n3-1;

a= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//ftiaxnoume to pinaka a pou 8a ginei o metasximatismos
A = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//ftiaxnoume to pinaka A pou einai o metasximatismenos
                                                            //exei xreisimopoii8ei to fftw_malloc tis fftw3.h ka8ws to xirizonte 
                                                            //kalitera oi sunartiseis tis fftw3.h
plan_a = fftw_plan_dft_2d(m,n, a, A, FFTW_FORWARD, FFTW_ESTIMATE);//ftiaxnoume to plano me orismata n0 n1 oi diastaseis
                                                                   //a h eisodos, A h eksodos kai FFTW_FORWARD o eu8us metasximatismos
                                                                  //etoimasame dld na ginei to a--F-->A, DEN exei ginei akoma 

b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//antistoixa kai edw
B = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);
plan_b = fftw_plan_dft_2d(m,n, b, B, FFTW_FORWARD, FFTW_ESTIMATE);

y= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//antistoixa kai edw me prosoxi omos Y eisodo kai y eksodo
Y = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//kai FFTW_BACKWARD gia na ektelestei o antistrofos
plan_y = fftw_plan_dft_2d(m,n, Y, y, FFTW_BACKWARD, FFTW_ESTIMATE);//Y--F(-1)-->y

insert(a,b,n0,n1,n2,n3,lim,n);//bazoume stoixia ston pinaka a kai ston b
printf("\nprinting a......\n");
if((f=fopen("rand_a","w"))==NULL){printf("ERROR OPENING FILE");exit(1);}
print(a,n0,n1,n,f);//ektiposi a
fclose(f);
printf("\nprinting b......\n");
if((f=fopen("rand_b","w"))==NULL){printf("ERROR OPENING FILE");exit(1);}
print(b,n2,n3,n,f);//ektiposi b
fclose(f);

zero(a,n0,n1,m,n);
zero(b,n2,n3,m,n);

choise='r';
// END EISODOS RAND complex
}
if (choise=='s')
{
// EISODOS RAND integer

printf("\n dose diastaseis n0xn1 n2xn3\n");
printf(" n0=");n0=input_dim();
printf(" n1=");n1=input_dim();
printf(" n2=");n2=input_dim();
printf(" n3=");n3=input_dim();
printf(" dose orio tuxaion ari8mon ");
lim=input_dim();

m=n0+n2-1;
n=n1+n3-1;

a= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//ftiaxnoume to pinaka a pou 8a ginei o metasximatismos
A = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//ftiaxnoume to pinaka A pou einai o metasximatismenos
                                                            //exei xreisimopoii8ei to fftw_malloc tis fftw3.h ka8ws to xirizonte 
                                                            //kalitera oi sunartiseis tis fftw3.h
plan_a = fftw_plan_dft_2d(m,n, a, A, FFTW_FORWARD, FFTW_ESTIMATE);//ftiaxnoume to plano me orismata n0 n1 oi diastaseis
                                                                   //a h eisodos, A h eksodos kai FFTW_FORWARD o eu8us metasximatismos
                                                                  //etoimasame dld na ginei to a--F-->A, DEN exei ginei akoma 

b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//antistoixa kai edw
B = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);
plan_b = fftw_plan_dft_2d(m,n, b, B, FFTW_FORWARD, FFTW_ESTIMATE);

y= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//antistoixa kai edw me prosoxi omos Y eisodo kai y eksodo
Y = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//kai FFTW_BACKWARD gia na ektelestei o antistrofos
plan_y = fftw_plan_dft_2d(m,n, Y, y, FFTW_BACKWARD, FFTW_ESTIMATE);//Y--F(-1)-->y

insert_int(a,b,n0,n1,n2,n3,lim,n);//bazoume stoixia ston pinaka a kai ston b
printf("\nprinting a......\n");
if((f=fopen("rand_a","w"))==NULL){printf("ERROR OPENING FILE");exit(1);}
print_real(a,n0,n1,n,f);//ektiposi a
fclose(f);
printf("\nprinting b......\n");
if((f=fopen("rand_b","w"))==NULL){printf("ERROR OPENING FILE");exit(1);}
print_real(b,n2,n3,n,f);//ektiposi b
fclose(f);


zero(a,n0,n1,m,n);
zero(b,n2,n3,m,n);


choise='s';
// END EISODOS RAND integer
}

if (choise=='f')
{
// EISODOS RAND real

printf("\n dose diastaseis n0xn1 n1xn2\n");
printf(" n0=");n0=input_dim();
printf(" n1=");n1=input_dim();
printf(" n2=");n2=input_dim();
printf(" n3=");n3=input_dim();
printf(" dose orio tuxaion ari8mon ");
lim=input_dim();

m=n0+n2-1;
n=n1+n3-1;

a= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//ftiaxnoume to pinaka a pou 8a ginei o metasximatismos
A = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//ftiaxnoume to pinaka A pou einai o metasximatismenos
                                                            //exei xreisimopoii8ei to fftw_malloc tis fftw3.h ka8ws to xirizonte 
                                                            //kalitera oi sunartiseis tis fftw3.h
plan_a = fftw_plan_dft_2d(m,n, a, A, FFTW_FORWARD, FFTW_ESTIMATE);//ftiaxnoume to plano me orismata n0 n1 oi diastaseis
                                                                   //a h eisodos, A h eksodos kai FFTW_FORWARD o eu8us metasximatismos
                                                                  //etoimasame dld na ginei to a--F-->A, DEN exei ginei akoma 

b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//antistoixa kai edw
B = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);
plan_b = fftw_plan_dft_2d(m,n, b, B, FFTW_FORWARD, FFTW_ESTIMATE);

y= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//antistoixa kai edw me prosoxi omos Y eisodo kai y eksodo
Y = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//kai FFTW_BACKWARD gia na ektelestei o antistrofos
plan_y = fftw_plan_dft_2d(m,n, Y, y, FFTW_BACKWARD, FFTW_ESTIMATE);//Y--F(-1)-->y

insert_real(a,b,n0,n1,n2,n3,lim,n);//bazoume stoixia ston pinaka a kai ston b
printf("\nprinting a......\n");
if((f=fopen("rand_a","w"))==NULL){printf("ERROR OPENING FILE");exit(1);}
print_real(a,n0,n1,n,f);//ektiposi a
fclose(f);
printf("\nprinting b......\n");
if((f=fopen("rand_b","w"))==NULL){printf("ERROR OPENING FILE");exit(1);}
print_real(b,n2,n3,n,f);//ektiposi b
fclose(f);


zero(a,n0,n1,m,n);
zero(b,n2,n3,m,n);


choise='s';
// END EISODOS RAND real
}

if (choise=='p')
// EISODOS FILE real
{
	
		if((f=fopen("input_a","r"))==NULL)                                // elenxos an uparxei to arxeio
	{printf("den brethike to arxeio 'input_a' 8a dimiourgi8ei kainourio\n");
	if((f=fopen("input_a","w"))==NULL){printf("ERROR OPENING FILE\n");exit(1);}
	fclose(f);
	if((f=fopen("input_a","r"))==NULL){printf("ERROR OPENING FILE\n");exit(1);}
	}
error=error_real(f); 
fclose(f);

while(error==1)                                                            //elenxos sta dedomena tou arxeiou
{printf(" la8os sta dedomena tou arxeiou 'input_a' dior8oste kai meta patiste enter....\n");
getchar();
if((f=fopen("input_a","r"))==NULL){printf("ERROR OPENING FILE\n");exit(1);}
error=error_real(f);
fclose(f);
}                            

		if((f=fopen("input_b","r"))==NULL)                                 // elenxos an uparxei to arxeio
	{printf("den brethike to arxeio 'input_b' 8a dimiourgi8ei kainourio\n");
	if((f=fopen("input_b","w"))==NULL){printf("ERROR OPENING FILE\n");exit(1);}
	fclose(f);
	if((f=fopen("input_b","r"))==NULL){printf("ERROR OPENING FILE\n");exit(1);}
	}
error=error_real(f); 
fclose(f);

while(error==1)                                                         //elenxos sta dedomena tou arxeiou
{printf(" la8os sta dedomena tou arxeiou 'input_b' dior8oste kai meta patiste enter....\n");
getchar();
if((f=fopen("input_b","r"))==NULL){printf("ERROR OPENING FILE\n");exit(1);}
error=error_real(f);
fclose(f);
}

if((f=fopen("input_a","r"))==NULL){printf("ERROR OPENING FILE");exit(1);}
lim=plithos(f);                             //briskoume to plithos twn stoixeiwn tou pinaka a
fclose(f);

if((f=fopen("input_a","r"))==NULL){printf("ERROR OPENING FILE");exit(1);}n0=0;
while(!feof(f))
{c = fgetc(f);
  if(c == '\n')
  {n0++;}                               //diabazei ka8e seira etsi mporoume na metrisoume tis seires
}
fclose(f);

n1=lim/n0;                          //etsi oi stoiles tou a einai oi pli8os/grammes

if((f=fopen("input_b","r"))==NULL){printf("ERROR OPENING FILE");exit(1);}
lim=plithos(f);                             //briskoume to plithos twn stoixeiwn tou pinaka b
fclose(f);

if((f=fopen("input_b","r"))==NULL){printf("ERROR OPENING FILE");exit(1);}n1_aux=0;
while(!feof(f))
{c = fgetc(f);
  if(c == '\n')
  {n1_aux++;}                  //diabazei ka8e seira etsi mporoume na metrisoume tis seires
}
fclose(f);

n2=lim/n1_aux;                                  //etsi oi stoiles tou b einai oi pli8os/grammes



printf("\n oi diastaseis einai %dx%d kai %dx%d\n",n0,n1,n1_aux,n2);

m=n0+n1_aux-1;
n=n1+n2-1;

printf("\n oi diastaseis tis suneliksis einai %dx%d\n",m,n);

a= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//ftiaxnoume to pinaka a pou 8a ginei o metasximatismos
A = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//ftiaxnoume to pinaka A pou einai o metasximatismenos
                                                            //exei xreisimopoii8ei to fftw_malloc tis fftw3.h ka8ws to xirizonte 
                                                            //kalitera oi sunartiseis tis fftw3.h
plan_a = fftw_plan_dft_2d(m,n, a, A, FFTW_FORWARD, FFTW_ESTIMATE);//ftiaxnoume to plano me orismata n0 n1 oi diastaseis
                                                                   //a h eisodos, A h eksodos kai FFTW_FORWARD o eu8us metasximatismos
                                                                  //etoimasame dld na ginei to a--F-->A, DEN exei ginei akoma 

b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//antistoixa kai edw
B = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);
plan_b = fftw_plan_dft_2d(m,n, b, B, FFTW_FORWARD, FFTW_ESTIMATE);

y= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//antistoixa kai edw me prosoxi omos Y eisodo kai y eksodo
Y = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//kai FFTW_BACKWARD gia na ektelestei o antistrofos
plan_y = fftw_plan_dft_2d(m,n, Y, y, FFTW_BACKWARD, FFTW_ESTIMATE);//Y--F(-1)-->y


if((f=fopen("input_a","r"))==NULL){printf("ERROR OPENING FILE");exit(1);}
insert_file(f,a,n0,n1,n);               //eisagogi a
fclose(f);

if((f=fopen("input_b","r"))==NULL){printf("ERROR OPENING FILE");exit(1);}
insert_file(f,b,n1_aux,n2,n);                  //isagogi b
fclose(f);
choise='p';

zero(a,n0,n1,m,n);
zero(b,n1_aux,n2,m,n);

// END EISODOS FILE real
}

if(choise=='c')
// EISODOS FILE complex
{
			if((f=fopen("input_a_complex","r"))==NULL)                                // elenxos an uparxei to arxeio
	{printf("den brethike to arxeio 'input_a_complex' 8a dimiourgi8ei kainourio\n");
	if((f=fopen("input_a_complex","w"))==NULL){printf("ERROR OPENING FILE\n");exit(1);}
	fclose(f);
	if((f=fopen("input_a_complex","r"))==NULL){printf("ERROR OPENING FILE\n");exit(1);}
	}
error=error_complex(f); 
fclose(f);

while(error==1)                                                            //elenxos sta dedomena tou arxeiou
{printf(" la8os sta dedomena tou arxeiou 'input_a_complex' dior8oste kai meta patiste enter....\n");
getchar();
if((f=fopen("input_a_complex","r"))==NULL){printf("ERROR OPENING FILE\n");exit(1);}
error=error_complex(f);
fclose(f);
}                            

		if((f=fopen("input_b_complex","r"))==NULL)                                 // elenxos an uparxei to arxeio
	{printf("den brethike to arxeio 'input_b_complex' 8a dimiourgi8ei kainourio\n");
	if((f=fopen("input_b_complex","w"))==NULL){printf("ERROR OPENING FILE\n");exit(1);}
	fclose(f);
	if((f=fopen("input_b_complex","r"))==NULL){printf("ERROR OPENING FILE\n");exit(1);}
	}
error=error_complex(f); 
fclose(f);

while(error==1)                                                        //elenxos sta dedomena tou arxeiou
{printf(" la8os sta dedomena tou arxeiou 'input_b_complex' dior8oste kai meta patiste enter....\n");
getchar();
if((f=fopen("input_b_complex","r"))==NULL){printf("ERROR OPENING FILE\n");exit(1);}
error=error_complex(f);
fclose(f);
}
	

if((f=fopen("input_a_complex","r"))==NULL){printf("ERROR OPENING FILE");exit(1);}
lim=plithos_complex(f);                             //briskoume to plithos twn stoixeiwn tou pinaka a
fclose(f);

if((f=fopen("input_a_complex","r"))==NULL){printf("ERROR OPENING FILE");exit(1);}n0=0;
while(!feof(f))
{c = fgetc(f);
  if(c == '\n')
  {n0++;}                             //diabazei ka8e seira etsi mporoume na metrisoume tis seires
}
fclose(f);

n1=lim/n0;                            //etsi oi stoiles tou a einai oi pli8os/grammes

if((f=fopen("input_b_complex","r"))==NULL){printf("ERROR OPENING FILE");exit(1);}
lim=plithos_complex(f);                             //briskoume to plithos twn stoixeiwn tou pinaka b
fclose(f);

if((f=fopen("input_b_complex","r"))==NULL){printf("ERROR OPENING FILE");exit(1);}n1_aux=0;
while(!feof(f))
{ c = fgetc(f);
  if(c == '\n')
  {n1_aux++;}                     //diabazei ka8e seira etsi mporoume na metrisoume tis seires

}
fclose(f);

n2=lim/n1_aux;                                 //etsi oi stoiles tou b einai oi pli8os/grammes



printf("\n oi diastaseis einai %dx%d kai %dx%d\n",n0,n1,n1_aux,n2);

m=n0+n1_aux-1;
n=n1+n2-1;

printf("\n oi diastaseis tis suneliksis einai %dx%d\n",m,n);

a= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//ftiaxnoume to pinaka a pou 8a ginei o metasximatismos
A = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//ftiaxnoume to pinaka A pou einai o metasximatismenos
                                                            //exei xreisimopoii8ei to fftw_malloc tis fftw3.h ka8ws to xirizonte 
                                                            //kalitera oi sunartiseis tis fftw3.h
plan_a = fftw_plan_dft_2d(m,n, a, A, FFTW_FORWARD, FFTW_ESTIMATE);//ftiaxnoume to plano me orismata n0 n1 oi diastaseis
                                                                   //a h eisodos, A h eksodos kai FFTW_FORWARD o eu8us metasximatismos
                                                                  //etoimasame dld na ginei to a--F-->A, DEN exei ginei akoma 

b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//antistoixa kai edw
B = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);
plan_b = fftw_plan_dft_2d(m,n, b, B, FFTW_FORWARD, FFTW_ESTIMATE);

y= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//antistoixa kai edw me prosoxi omos Y eisodo kai y eksodo
Y = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*m*n);//kai FFTW_BACKWARD gia na ektelestei o antistrofos
plan_y = fftw_plan_dft_2d(m,n, Y, y, FFTW_BACKWARD, FFTW_ESTIMATE);//Y--F(-1)-->y


if((f=fopen("input_a_complex","r"))==NULL){printf("ERROR OPENING FILE");exit(1);}
insert_file_complex(f,a,n0,n1,n);               //eisagogi a
fclose(f);

if((f=fopen("input_b_complex","r"))==NULL){printf("ERROR OPENING FILE");exit(1);}
insert_file_complex(f,b,n1_aux,n2,n);                  //eisagogi b
fclose(f);

zero(a,n0,n1,m,n);
zero(b,n1_aux,n2,m,n);

choise='c';

// END EISODOS FILE complex
}

clock_t start=clock();
//fftw_execute(plan_a);// autto 8a ektelousame xwris threads gia na broume to a--F-->A
//fftw_execute(plan_b); //kai b--F-->B antistoixa

thread_data_array[0].plan=plan_a;// to orisma tou protou thread einai to plano tou a
thread_data_array[1].plan=plan_b;//tou deuterou to plano tou b
rc=pthread_create(&threads[0],&attr,execute,(void *) &thread_data_array[0]);//ftiaxnoume ena thread pou ektelei to a--F-->A
rc=pthread_create(&threads[1],&attr,execute,(void *) &thread_data_array[1]);//kai ena pou ektelei to b--F-->B
pthread_attr_destroy(&attr);
rc=pthread_join(threads[0],&status);
rc=pthread_join(threads[1],&status);// telionoun ta thread


mult(Y,A,B,m,n);//ypologizoume to Y pou einai o polaplasiasmos simion A*B
fftw_execute(plan_y);//ekteloume to Y--F(-1)-->y
normalize_2D(y,m,n);//diairoume ta panta me m*n
clock_t end1=clock();float seconds=(float)(end1-start)/CLOCKS_PER_SEC;
printf("\nprinting y.......\n");




if((choise=='r')||(choise=='c'))
{
if((f=fopen("output","w"))==NULL){printf("ERROR OPENING FILE");exit(1);}
print(y,m,n,n,f);// kai ektuponoume to y pou einai pleon i suneliksi tou a me to b
fclose(f);
}

if((choise=='s')||(choise=='p'))
{
if((f=fopen("output","w"))==NULL){printf("ERROR OPENING FILE");exit(1);}
print_real(y,m,n,n,f);// kai ektuponoume to y pou einai pleon i suneliksi tou a me to b
fclose(f);
}

fftw_destroy_plan(plan_a);
fftw_destroy_plan(plan_b);
fftw_destroy_plan(plan_y);
fftw_free(a);fftw_free(A);
fftw_free(b);fftw_free(B);
fftw_free(y);fftw_free(Y);//sbinoume oti exoume dimiourgisei me tis sunartuseis tis fftw3.h
printf("xronos ektelesis suneliksis %f \n",seconds);choise='*';

}
}
pthread_exit(NULL);// eksodos asfaleias twn thread
}


