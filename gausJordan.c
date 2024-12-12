#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_DIM 15
#define EPSILON 1e-10

// Variáveis globais para a matriz aumentada e ordem do sistema
double matriz[MAX_DIM][MAX_DIM + 1];
int ordem;

// Função para encontrar o pivô de maior valor absoluto
void encontraPivoCompleto(int k, int *linha_pivo, int *col_pivo) {
    double max = 0.0;
    *linha_pivo = k;
    *col_pivo = k;

    for (int i = k; i < ordem; i++) {
        for (int j = k; j < ordem; j++) {
            if (fabs(matriz[i][j]) > max) {
                max = fabs(matriz[i][j]);
                *linha_pivo = i;
                *col_pivo = j;
            }
        }
    }
}

// Função para trocar duas linhas da matriz
void trocaLinhas(int linha1, int linha2) {
    for (int j = 0; j <= ordem; j++) {
        double temp = matriz[linha1][j];
        matriz[linha1][j] = matriz[linha2][j];
        matriz[linha2][j] = temp;
    }
}

// Função para trocar duas colunas da matriz
void trocaColunas(int col1, int col2) {
    for (int i = 0; i < ordem; i++) {
        double temp = matriz[i][col1];
        matriz[i][col1] = matriz[i][col2];
        matriz[i][col2] = temp;
    }
}

// Função para imprimir a matriz aumentada
void imprimeMatriz() {
    printf("\nMatriz aumentada:\n");
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j <= ordem; j++) {
            printf("%10.4f ", matriz[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// Função principal do método de Gauss-Jordan
int gaussJordan() {
    int *permutacao = (int *)malloc(ordem * sizeof(int));

    // Inicializa vetor de permutação
    for (int i = 0; i < ordem; i++) {
        permutacao[i] = i;
    }

    // Para cada coluna k
    for (int k = 0; k < ordem; k++) {
        // Encontra o pivô
        int linha_pivo, col_pivo;
        encontraPivoCompleto(k, &linha_pivo, &col_pivo);

        // Verifica se o pivô é nulo
        if (fabs(matriz[linha_pivo][col_pivo]) < EPSILON) {
            free(permutacao);
            return 0; // Sistema sem solução única
        }

        // Troca linhas e colunas se necessário
        if (linha_pivo != k) {
            trocaLinhas(k, linha_pivo);
        }
        if (col_pivo != k) {
            trocaColunas(k, col_pivo);
            int temp = permutacao[k];
            permutacao[k] = permutacao[col_pivo];
            permutacao[col_pivo] = temp;
        }

        // Normaliza a linha do pivô
        double pivo = matriz[k][k];
        for (int j = k; j <= ordem; j++) {
            matriz[k][j] /= pivo;
        }

        // Zera os elementos acima e abaixo do pivô
        for (int i = 0; i < ordem; i++) {
            if (i != k) {
                double fator = matriz[i][k];
                for (int j = k; j <= ordem; j++) {
                    matriz[i][j] -= fator * matriz[k][j];
                }
            }
        }

        printf("\nIteração %d:\n", k + 1);
        imprimeMatriz();
    }

    // Reordena a solução de acordo com as permutações realizadas
    double *solucao_temp = (double *)malloc(ordem * sizeof(double));
    for (int i = 0; i < ordem; i++) {
        solucao_temp[i] = matriz[i][ordem];
    }

    for (int i = 0; i < ordem; i++) {
        matriz[permutacao[i]][ordem] = solucao_temp[i];
    }

    free(solucao_temp);
    free(permutacao);
    return 1; // Sistema resolvido com sucesso
}

// Função para ler o sistema linear de um arquivo
int leSistema(const char *nomeArquivo) {
    FILE *arquivo = fopen(nomeArquivo, "r");
    if (!arquivo) {
        printf("Erro ao abrir o arquivo!\n");
        return 0;
    }

    // Lê a ordem do sistema
    fscanf(arquivo, "%d", &ordem);

    if (ordem > MAX_DIM) {
        printf("Ordem do sistema maior que o permitido!\n");
        fclose(arquivo);
        return 0;
    }

    // Lê os coeficientes da matriz
    for (int i = 0; i < ordem; i++) {
        for (int j = 0; j < ordem; j++) {
            fscanf(arquivo, "%lf", &matriz[i][j]);
        }
    }

    // Lê o vetor de resultados
    for (int i = 0; i < ordem; i++) {
        fscanf(arquivo, "%lf", &matriz[i][ordem]);
    }

    fclose(arquivo);
    return 1;
}

int main() {
    // Lê o sistema do arquivo
    if (!leSistema("Sistema.txt")) {
        return 1;
    }

    printf("Sistema original:\n");
    imprimeMatriz();

    // Resolve o sistema
    int resultado = gaussJordan();

    if (resultado) {
        printf("\nSolução do sistema:\n");
        for (int i = 0; i < ordem; i++) {
            printf("x%d = %.8f\n", i + 1, matriz[i][ordem]);
        }
    } else {
        printf("\nO sistema não possui solução única!\n");
    }

    return 0;
}
