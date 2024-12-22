#include <stdio.h>
#include <stdlib.h>
#include <math.h> // Necessário para fabs
#include<stdbool.h>

#define MAX_DIM 15
// Variáveis globais para a matriz aumentada e ordem do sistema
    long double matriz[MAX_DIM][MAX_DIM + 1]={0}; // Matriz aumentada (inclui a coluna do vetor de resultados)

int ordem=0; // Ordem do sistema (número de equações e variáveis)
int max_iter=0; // Máximo de iterações
long double epsilon=0; // Precisão para o critério de parada
int posicoes_linhas[MAX_DIM]={0}; // Armazena a posição original das linhas após as trocas de linhas
int posicoes_colunas[MAX_DIM]={0}; // Armazena a posição original das colunas após as trocas de colunas
//* OK
int verificaCriterioSassenfeld() {
    long double beta[MAX_DIM] = {0.0};
    for (int i = 0; i < ordem; i++) {
        long double soma_beta = 0.0;
        for (int j = 0; j < ordem; j++) {
            if (i != j) {
                // Para cada linha, calcula a soma dos coeficientes 
                // divididos pelo coeficiente diagonal, considerando os betas anteriores
                if (j < i) {
                    soma_beta += (fabsl(matriz[i][j]) / fabsl(matriz[i][i])) * beta[j];
                } else {
                    soma_beta += fabsl(matriz[i][j]) / fabsl(matriz[i][i]);
                }
            }
        }
        beta[i] = soma_beta;
        if (beta[i] >= 1.0) {
            printf("Critério de Sassenfeld: Método pode NÃO convergir!\n");
            return 0;
        }
    }
    printf("Critério de Sassenfeld: Método tem GARANTIA de convergência!\n");
    return 1;
}

int leSistema(const char *nomeArquivo) {
    FILE *arquivo = fopen(nomeArquivo, "r");
    if (!arquivo) {
        printf("Erro ao abrir o arquivo!\n");
        return 0;
    }
    fscanf(arquivo, "%d", &ordem);
    if (ordem > MAX_DIM) {
        printf("Ordem do sistema maior que o permitido!\n");
        fclose(arquivo);
        return 0;
    }
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j < ordem; j++) {
            fscanf(arquivo, "%Lf", &matriz[i][j]);
        }
    }
    for (int i = 0; i < ordem; i++) {
        fscanf(arquivo, "%Lf", &matriz[i][ordem]);
    }
    fscanf(arquivo, "%d %Lf", &max_iter, &epsilon);
    fclose(arquivo);
    for (int i = 0; i < ordem; i++) {
        posicoes_linhas[i] = i; // Inicialmente as linhas estão na posição original
        posicoes_colunas[i] = i; // Inicialmente as colunas estão na posição original
    }
    return 1;
}

//* OK
void imprimir_matriz() {
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j <= ordem; j++) {
            printf("%8.6Lf ", matriz[i][j]);
        }
        printf("\n");
    }
}
//* OK 
void encontraPivoCompleto(int k, int *linha_pivo, int *col_pivo) {
    long double max = 0.0;
    *linha_pivo = k;
    *col_pivo = k;

    // Percorre TODA a matriz, não apenas a submatriz
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j < ordem; j++) {
            if (fabsl(matriz[i][j]) > max) {
                max = fabsl(matriz[i][j]);
                *linha_pivo = i;
                *col_pivo = j;
            }
        }
    }
}
//* OK
void trocaLinhas(int linha1, int linha2) {
    for (int j = 0; j <= ordem; j++) {
        long double temp = matriz[linha1][j];
        matriz[linha1][j] = matriz[linha2][j];
        matriz[linha2][j] = temp;
    }
    
    int temp_pos = posicoes_linhas[linha1];
    posicoes_linhas[linha1] = posicoes_linhas[linha2];
    posicoes_linhas[linha2] = temp_pos;
}
//* OK
void trocaColunas(int col1, int col2) {
    for (int i = 0; i < ordem; i++) {
        long double temp = matriz[i][col1];
        matriz[i][col1] = matriz[i][col2];
        matriz[i][col2] = temp;
    }
    
    int temp_pos = posicoes_colunas[col1];
    posicoes_colunas[col1] = posicoes_colunas[col2];
    posicoes_colunas[col2] = temp_pos;
}


void gauss_seidel() {
    long double atual[MAX_DIM] = {0};
    long double anterior[MAX_DIM] = {0};
    int interacoes = 0;
    long double erro_max_atual = 0.0;
    // Aplicar pivotamento e preparação da matriz
    for (int k = 0; k < ordem; k++) {
        int linha_pivo, col_pivo;
        encontraPivoCompleto(k, &linha_pivo, &col_pivo);
        if (linha_pivo != k) trocaLinhas(k, linha_pivo);
        if (col_pivo != k) trocaColunas(k, col_pivo);
    }
  
    imprimir_matriz();
    while (interacoes < max_iter) {
        printf("\n--- Iteração %d ---\n", interacoes + 1);
        
        // Salva a solução anterior
        for (int i = 0; i < ordem; i++) {
            anterior[i] = atual[i];
        }

        // Método de Gauss-Seidel
        for (int i = 0; i < ordem; i++) {
            long double soma = 0.0;
            

            // Cálculo da solução
            for (int j = 0; j < ordem; j++) {
                if (i != j) {
                    soma += matriz[i][j] * atual[j];
                }
            }

            if (fabsl(matriz[i][i]) > 1e-15) {
                atual[i] = (matriz[i][ordem] - soma) / matriz[i][i];
                printf("x%d = (%.4Lf - %.4Lf) / %.4Lf = %.10Lf\n", 
                       i + 1, matriz[i][ordem], soma, matriz[i][i], atual[i]);
            } else {
                printf("\nErro: divisão por zero na posição [%d][%d] da matriz.\n", i, i);
                return;
            }
        }

        // Calcula o critério de erro máximo
        erro_max_atual = 0.0;
        printf("\nAnálise de erro:\n");
        for (int i = 0; i < ordem; i++) {
            long double erro_relativo = fabsl(atual[i] - anterior[i]);
            printf("  Erro relativo de x%d: %.10Lf\n", i + 1, erro_relativo);
            if (erro_relativo > erro_max_atual) {
                erro_max_atual = erro_relativo;
            }
        }

        printf("  Erro máximo da iteração %d: %.10Lf\n", interacoes + 1, erro_max_atual);
        
        // Critério de parada por erro
        if (erro_max_atual < epsilon) {
            printf("\nConvergência atingida após %d iterações com erro máximo de %.10Lf.\n", 
                   interacoes + 1, erro_max_atual);
            break;
        }

        interacoes++;
    }

    // Verificações finais
    if (interacoes == max_iter) {
        printf("\nNúmero máximo de iterações atingido sem convergência.\n");
    }

    // Imprime a solução final
    printf("\nSolução final:\n");
    for (int i = 0; i < ordem; i++) {
        printf("x%d = %.10Lf (Posição Original: %d)\n", 
               i + 1, atual[i], posicoes_colunas[i] + 1);
    }
}
void escalonarMatriz() {
    for (int k = 0; k < ordem; k++) {
        // Encontra o maior elemento na coluna k (pivotamento parcial)
        int maxIndex = k;
        for (int i = k + 1; i < ordem; i++) {
            if (fabsl(matriz[i][k]) > fabsl(matriz[maxIndex][k])) {
                maxIndex = i;
            }
        }

        // Troca a linha k com a linha maxIndex
        if (maxIndex != k) {
            trocaLinhas(k, maxIndex);
        }

        // Elimina os elementos abaixo do pivô
        for (int i = k + 1; i < ordem; i++) {
            if (fabsl(matriz[k][k]) > epsilon) {
                long double fator = matriz[i][k] / matriz[k][k];
                for (int j = k; j <= ordem; j++) { // Inclui o termo independente na matriz aumentada
                    matriz[i][j] -= fator * matriz[k][j];
                }
            }
        }
    }
}

void classificarSistema() {
    escalonarMatriz(); // Escalona a matriz antes de classificar
    int linhasNulas = 0;
    int rank = 0;

    for (int i = 0; i < ordem; i++) {
        bool linhaNulaCoeficientes = true;
        bool linhaNulaCompleta = true;

        for (int j = 0; j < ordem; j++) {
            if (fabsl(matriz[i][j]) > epsilon) {
                linhaNulaCoeficientes = false;
                break;
            }
        }

        for (int j = 0; j <= ordem; j++) { // Inclui o termo independente na verificação
            if (fabsl(matriz[i][j]) > epsilon) {
                linhaNulaCompleta = false;
                break;
            }
        }

        if (linhaNulaCoeficientes && !linhaNulaCompleta) {
            // Linha de coeficientes nula, mas termo independente não nulo
            printf("O sistema é classificado como SI (Sistema Impossível).\n");
            return;
        } 

        if (!linhaNulaCoeficientes) {
            rank++;
        }
    }

    if (rank < ordem) {
        if (rank < ordem - 1) {
            printf("O sistema é classificado como SPI (Sistema Possível e Indeterminado).\n");
        } else {
            printf("O sistema é classificado como SPI (Sistema Possível e Indeterminado).\n");
        }
    } else {
        printf("O sistema é classificado como SPD (Sistema Possível e Determinado).\n");
    }
}
int main() {
    if (leSistema("../Ex1.txt")) {
        printf("Sistema lido com sucesso!\n");
        escalonarMatriz();
        if (verificaCriterioSassenfeld()) {
            gauss_seidel(); // Chama o método de Gauss-Seidel apenas se o critério for satisfeito
        } else {
            printf("Não é possível garantir a convergência do método de Gauss-Seidel.\n");
            gauss_seidel();
        }// Chama o método de Gauss-Seidel
        classificarSistema();
    } else {
        printf("Falha ao ler o sistema.\n");
    }
    return 0;
}
