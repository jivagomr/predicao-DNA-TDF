# predicao-DNA-TDF
Predição de zonas codificadores em cadeias de DNA utilizando Transformada Discreta de Fourier.

Alogoritmo desenvolvido para o trabalho final da disciplina de Biologia Computacional ministrada pela profa. Helena Cristina da Gama Leitão no programa de pós-graduação da UFF em 2011. 

Fickett 1982 apresentou um estudo onde são listadas diversas características de regiões codificadoras e não codificadoras. Uma dessas característica forneceu a base necessária para a utilização de processamento de sinais digitais como ferramenta na predição de genes; a
periodicidade três das regiões codificadoras. Nas regiões codificadoras, uma determinada base se repete com período 3 (ou frequência = 1/3), característica essa que não está presente em regiões não codificadoras. A presença de periodicidade 3 em regiões codificadoras e a ausência em regiões não codificadoras é considerada uma característica universal das sequências de DNA (Tiware et al. 1997).

Voss 1992 demonstrou que era possível observar a característica da periodicidade 3 utilizando Análise de Fourier em uma sequência de DNA. A característica estaria visível por meio de um pico espectral na frequência 1/3.

A implementação foi feita utilizando a linguagem C e compilada utilizando o gcc (http://gcc.gnu.org/) em um ambiente Linux utilizando a distribuição Fedora 15 (http://fedoraproject.org/) com arquitetura x86_64 e kernel do Linux na versão 2.6.38.8. As Transformadas de Fourier foram feitas utilizando a API FFTW na versão 3.2.2 (Fastest Fourier Transform in the West), também implementada em C e tida como uma das das implementações mais rápidas de FFT. Para compilar o programa é necessário “linkar” a biblioteca relativa a FFTW, como por exemplo: gcc -o fft-jivago -lfftw3 -lm fft-jivago.c Após compilado, para executar o programa é necessário passar por parâmetro o arquivo no padrão GenBank contendo a sequência a ser analisada. Como por exemplo: ./fft-jivago sequences/AF135270.gb.

REFERÊNCIAS 

Fickett. J., W. Recognition of protein coding regions in DNA sequences. Nucleic Acids Res., 1982. 10, 5303-5318.

Pessoa, S., L. Análise da informação mútua em sequências de DNA homólogas. Dissertação de Mestrado – Programa de Pós Graduação em Computação, Universidade Federal Fluminense, Niterói – RJ. 2004.

Tiwari S., Ramachandran S., Bhattacharya A., Bhattacharya S., Ramaswamy R. Prediction of probable genes by Fourier analysis of genomic sequences. Comput. Appl. Biosci. 1997; 13 (3) : 263-270

Voss., R., F.Evolution of long-range fractal correlations and l/f noise in DNA base sequences. Phys. Rev. Lett., 1992. 68, 3805-3808
