#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <math.h>

struct base {
    int nct;
    struct base *prox;
};

struct base *inicio_SEQ;
struct base *fim_SEQ;
int tamanho;

void inicia_seq() {
    inicio_SEQ=NULL;
    fim_SEQ=NULL;
}

int add_base(int nct){

    struct base *novo_nct;
    //struct base *aux;

    if (( novo_nct = (struct base *) malloc(sizeof(struct base))) == NULL) {
        puts( "Falta Memoria\n");  return 0;
    }
    novo_nct->nct = nct;
    novo_nct->prox = NULL;
 
    if (inicio_SEQ == NULL) {
        inicio_SEQ = novo_nct;
	fim_SEQ = inicio_SEQ;
    }
    else{
	fim_SEQ->prox=novo_nct;
	fim_SEQ=novo_nct;
    }
    return 1;
}

int expected_origin(int cur) {

	int expected=-1;
	
	//ORIGIN
        //0-> 79 	O
	//1-> 82 	R
	//2-> 73 	I
	//3-> 71 	G
	//4-> 73 	I
	//5-> 78 	N

	switch (cur) {
		case 0:
			expected=79;
			break;
		case 1:
			expected=82;
			break;
		case 2:
			expected=73;
			break;
		case 3:
			expected=71;
			break;
		case 4:
			expected=73;
			break;
		case 5:
			expected=78;
			break;
	}


	return expected;

}

void ler_sequencia(char *seq) {
	
	//a 97
	//c 99
	//g 103
	//t 116

	int origin=0;
	int expected;
	int seq_len=0;

	FILE *fp;
	fp=fopen(seq, "r");

	char ch;
	int ch_code;

	int t;

	if (fp != NULL) {
		inicia_seq(&inicio_SEQ); //inicia a sequencia de DNA (lista dinâmica)
		ch = getc(fp);
		ch_code = ch;
		while (ch!=EOF)  {
			ch = getc(fp);
			ch_code = ch;
			if (origin==6) {
				if ((ch_code==97) || (ch_code==99) || (ch_code==103) || (ch_code==116)) {
					//printf("%c",ch);
					add_base(ch_code);
					seq_len++;
				}
			}
			else {
				expected = expected_origin(origin);
				if (ch_code==expected)
					origin++;
				else
					origin=0;
			}
        	}
		tamanho=seq_len;
		printf("\nTamanho da sequencia: %d\n",seq_len);
		double f;
		f=ceil(tamanho/3);
		printf("\nPico esperado em k=N/3, com K = %.2f\n",f);
	}
	else {
		printf("nao foi possivel ler arquivo");
	}


}

double calcKvalue(double num) {
	double calc;	
	calc = (num<0) ? (num*(-1)) : num;
	calc = (num*num);
	return calc;
}

void executar(char *out_file) {
	int Ua[tamanho];
	int Uc[tamanho];
	int Ug[tamanho];
	int Ut[tamanho];

	struct base *atual;
	struct base *aux;

	int nct;

	int x,i,j;

	atual = inicio_SEQ;

	x=0;

	while(atual != NULL){
		nct = atual->nct;
	
		if (nct==97) {
			Ua[x]=1;
			Uc[x]=0;
			Ug[x]=0;
			Ut[x]=0;
		}
		else if (nct==99) {
			Ua[x]=0;
			Uc[x]=1;
			Ug[x]=0;
			Ut[x]=0;
		}
		else if (nct==103) {
			Ua[x]=0;
			Uc[x]=0;
			Ug[x]=1;
			Ut[x]=0;
		}
		else if (nct==116) {
			Ua[x]=0;
			Uc[x]=0;
			Ug[x]=0;
			Ut[x]=1;
		}

		aux = atual;	    
		atual = atual->prox;
		free(aux);
		x++;
	}

	//AGORA EU FAÇO A TRANSFORMADA DOS 	
	fftw_complex *in, *out,*fftUa,*fftUc,*fftUg,*fftUt,*fftSequence;
	fftw_plan p;

	fftUa = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tamanho);
	fftUc = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tamanho);
	fftUg = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tamanho);
	fftUt = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tamanho);

	fftSequence = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tamanho);

	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tamanho);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * tamanho);
	p = fftw_plan_dft_1d(tamanho, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	//FFT Ua
	// passado para in o valor do vetor Ua
	for (i = 0; i < tamanho; i++) {
		in[i][0] = Ua[i];
		in[i][1] = 0;
	}
	fftw_execute(p); /* executando FFT para Ua */
	
	//devolvendo pra fftUa o valor de out
	for (i = 0; i < tamanho; i++) {
		fftUa[i][0] = out[i][0];
		fftUa[i][1] = out[i][1];
	}

	//FFT Uc
	// passado para in o valor do vetor Ua
	for (i = 0; i < tamanho; i++) {
		in[i][0] = Uc[i];
		in[i][1] = 0;
	}
    fftw_execute(p); /* executando FFT para Ua */
	   
	//devolvendo pra fftUa o valor de out
	for (i = 0; i < tamanho; i++) {
		fftUc[i][0] = out[i][0];
		fftUc[i][1] = out[i][1];
	}

	//FFT Ug
	// passado para in o valor do vetor Ua
	for (i = 0; i < tamanho; i++) {
		in[i][0] = Ug[i];
		in[i][1] = 0;
	}
    fftw_execute(p); /* executando FFT para Ua */
	
	//devolvendo pra fftUa o valor de out
	for (i = 0; i < tamanho; i++) {
		fftUg[i][0] = out[i][0];
		fftUg[i][1] = out[i][1];
	}

	//FFT Ut
	// passado para in o valor do vetor Ua
	for (i = 0; i < tamanho; i++) {
		in[i][0] = Ut[i];
		in[i][1] = 0;
	}
	fftw_execute(p); /* executando FFT para Ua */
	
	//devolvendo pra fftUa o valor de out
	for (i = 0; i < tamanho; i++) {
		fftUt[i][0] = out[i][0];
		fftUt[i][1] = out[i][1];
	}

	int int_part;


	//FAZENDO O SOMATÓRIO DAS FFTs
	double kA,kC,kG,kT;
	int n;

	n=ceil(tamanho/2);

	for (i = 1; i < n; i++) {

		//parte real
		kA=calcKvalue(fftUa[i][0]);
		kC=calcKvalue(fftUc[i][0]);
		kG=calcKvalue(fftUg[i][0]);
		kT=calcKvalue(fftUt[i][0]);
		fftSequence[i][0]=kA+kC+kG+kT;

		//parte imaginária		
		kA=calcKvalue(fftUa[i][1]);
		kC=calcKvalue(fftUc[i][1]);
		kG=calcKvalue(fftUg[i][1]);
		kT=calcKvalue(fftUt[i][1]);
		fftSequence[i][1]=kA+kC+kG+kT;

	}

	//ESCREVENDO A SAÍDA NO ARQUIVO
	FILE *fp_saida;

	fp_saida = fopen ("processadas/fft-sequence.txt", "w+");
	if (fp_saida != NULL) {
		for (i = 1; i < n; i++) {
			//fputs (realPart, fp_saida);
			//fprintf(fp_saida,"%f %f \n",fftSequence[i][0],fftSequence[i][1]);
			fprintf(fp_saida,"%i %f \n",(i+1),fftSequence[i][0]);
			//printf("%f",out[i][0]);
		}
		fclose(fp_saida);
	}

	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	fftw_free(fftUa);
	fftw_free(fftUc);
	fftw_free(fftUg);
	fftw_free(fftUt);

}

void gnuplot(const char *gnucommand) {
  char syscommand[1024];
  sprintf(syscommand, "echo \"%s\" | gnuplot -persist", gnucommand);
  system(syscommand);
}

int main (int argc, char *argv[]) {

	char nct;
	
	ler_sequencia(argv[1]);
	executar(argv[1]);

	gnuplot("plot 'processadas/fft-sequence.txt' with lines");

	return 0;
}