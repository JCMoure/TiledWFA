#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX_LONG 100
#define MATCH 0
#define MISMATCH 1
#define GAP 1

// Estructura para almacenar el resultado del alineamiento
typedef struct {
    char seq1_alineada[MAX_LONG * 2];
    char seq2_alineada[MAX_LONG * 2];
    int puntuacion;
} Alineamiento;

int minimo(int a, int b, int c) {
    int min = a;
    if (b < min) min = b;
    if (c < min) min = c;
    return min;
}

int obtener_puntuacion(char a, char b) {    
    if (a == b) {
        return MATCH;
    } else {
        return MISMATCH;
    }
}

// Versión modificada de alinear_secuencias que puede mostrar las matrices (Needleman-Wunsch)
Alineamiento alinear_secuencias_con_matriz(char *seq1, char *seq2, int mostrar_matriz) {
    int len1 = strlen(seq1);
    int len2 = strlen(seq2);
    
    // Matriz de puntuaciones
    int matriz[MAX_LONG][MAX_LONG];
    // Matriz para trazado (0: diagonal, 1: arriba, 2: izquierda)
    int traza[MAX_LONG][MAX_LONG];
    
    // Inicializar primera fila y primera columna
    matriz[0][0] = 0;
    traza[0][0] = 0;
    
    for (int i = 1; i <= len1; i++) {
        matriz[i][0] = matriz[i-1][0] + GAP;
        traza[i][0] = 1; // Viene de arriba
    }
    
    for (int j = 1; j <= len2; j++) {
        matriz[0][j] = matriz[0][j-1] + GAP;
        traza[0][j] = 2; // Viene de izquierda
    }
    
    // Llenar la matriz
    for (int i = 1; i <= len1; i++) {
        for (int j = 1; j <= len2; j++) {
            int match = matriz[i-1][j-1] + obtener_puntuacion(seq1[i-1], seq2[j-1]);
            int del = matriz[i-1][j] + GAP;     // Gap en seq2
            int ins = matriz[i][j-1] + GAP;     // Gap en seq1
            
            matriz[i][j] = minimo(match, del, ins);
            
            // Guardar dirección para trazado
            if (matriz[i][j] == match) {
                traza[i][j] = 0; // Diagonal
            } else if (matriz[i][j] == del) {
                traza[i][j] = 1; // Arriba
            } else {
                traza[i][j] = 2; // Izquierda
            }
        }
    }
    
    // Mostrar matrices si se solicita
    if (mostrar_matriz) {
        mostrar_matriz_con_traza(matriz, traza, len1, len2, seq1, seq2);
    }
    
    Alineamiento resultado;
    resultado.puntuacion = matriz[len1][len2];
    
    int i = len1, j = len2;
    int pos1 = 0, pos2 = 0;
    
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && traza[i][j] == 0) {
            resultado.seq1_alineada[pos1++] = seq1[i-1];
            resultado.seq2_alineada[pos2++] = seq2[j-1];
            i--; j--;
        } else if (i > 0 && traza[i][j] == 1) {
            resultado.seq1_alineada[pos1++] = seq1[i-1];
            resultado.seq2_alineada[pos2++] = '-';
            i--;
        } else if (j > 0 && traza[i][j] == 2) {
            resultado.seq1_alineada[pos1++] = '-';
            resultado.seq2_alineada[pos2++] = seq2[j-1];
            j--;
        }
    }
    
    // Invertir las secuencias
    resultado.seq1_alineada[pos1] = '\0';
    resultado.seq2_alineada[pos2] = '\0';
    
    // Invertir las cadenas
    for (int k = 0; k < pos1/2; k++) {
        char temp = resultado.seq1_alineada[k];
        resultado.seq1_alineada[k] = resultado.seq1_alineada[pos1-1-k];
        resultado.seq1_alineada[pos1-1-k] = temp;
    }
    
    for (int k = 0; k < pos2/2; k++) {
        char temp = resultado.seq2_alineada[k];
        resultado.seq2_alineada[k] = resultado.seq2_alineada[pos2-1-k];
        resultado.seq2_alineada[pos2-1-k] = temp;
    }
    
    return resultado;
}

void mostrar_matriz_con_traza(int matriz[MAX_LONG][MAX_LONG], 
                              int traza[MAX_LONG][MAX_LONG], 
                              int len1, int len2, char *seq1, char *seq2) {
    printf("\n=== MATRIZ DE COSTES CON TRAZADO ===\n");
    printf("Filas: Secuencia 1 (%s)\n", seq1);
    printf("Columnas: Secuencia 2 (%s)\n", seq2);
    printf("Formato: [coste|dirección] (←:izquierda, ↑:arriba, ↖:diagonal)\n\n");
    
    // Imprimir encabezado de columnas
    printf("        ");  // Espacio para la esquina
    printf("   -    ");  // Gap inicial
    for (int j = 0; j < len2; j++) {
        printf("   %c    ", seq2[j]);
    }
    printf("\n");
    
    // Imprimir línea separadora
    printf("  ");
    for (int j = 0; j <= len2 + 1; j++) {
        printf("--------");
    }
    printf("\n");
    
    // Imprimir cada fila
    for (int i = 0; i <= len1; i++) {
        // Encabezado de fila
        if (i == 0) {
            printf("- |");  // Gap inicial
        } else {
            printf("%c |", seq1[i-1]);
        }
        
        // Imprimir valores con dirección
        for (int j = 0; j <= len2; j++) {
            char direccion;
            if (i == 0 && j == 0) {
                direccion = '*';  // Punto de inicio
            } else if (traza[i][j] == 0) {
                direccion = '↖';  // Diagonal
            } else if (traza[i][j] == 1) {
                direccion = '↑';  // Arriba
            } else {
                direccion = '←';  // Izquierda
            }
            printf(" %3d %c ", matriz[i][j], direccion);
        }
        printf("\n");
        
        // Línea separadora
        printf("  |");
        for (int j = 0; j <= len2; j++) {
            printf("--------");
        }
        printf("\n");
    }
    printf("\n");
}

// Función para mostrar el alineamiento de forma visual
void mostrar_alineamiento(Alineamiento alineamiento) {
    int len = strlen(alineamiento.seq1_alineada);
    
    printf("\n=== RESULTADO DEL ALINEAMIENTO ===\n");
    printf("Puntuación total: %d\n\n", alineamiento.puntuacion);
    
    // Mostrar secuencia 1
    printf("Seq1: ");
    for (int i = 0; i < len; i++) {
        printf("%c ", alineamiento.seq1_alineada[i]);
    }
    printf("\n");
    
    // Mostrar indicadores de match/mismatch
    printf("      ");
    for (int i = 0; i < len; i++) {
        if (alineamiento.seq1_alineada[i] == alineamiento.seq2_alineada[i]) {
            printf("| ");  // Match
        } else if (alineamiento.seq1_alineada[i] == '-' || 
                   alineamiento.seq2_alineada[i] == '-') {
            printf("  ");  // Gap
        } else {
            printf("* ");  // Mismatch
        }
    }
    printf("\n");
    
    // Mostrar secuencia 2
    printf("Seq2: ");
    for (int i = 0; i < len; i++) {
        printf("%c ", alineamiento.seq2_alineada[i]);
    }
    printf("\n");
    
    // Calcular y mostrar porcentaje de identidad
    int matches = 0;
    int total = 0;
    for (int i = 0; i < len; i++) {
        if (alineamiento.seq1_alineada[i] != '-' && 
            alineamiento.seq2_alineada[i] != '-') {
            total++;
            if (alineamiento.seq1_alineada[i] == alineamiento.seq2_alineada[i]) {
                matches++;
            }
        }
    }
    
    if (total > 0) {
        float identidad = (float)matches / total * 100;
        printf("\nIdentidad: %.1f%% (%d/%d posiciones coincidentes)\n", 
               identidad, matches, total);
    }
}

int main(int argc, char **argv) {
    char secuencia1[MAX_LONG];
    char secuencia2[MAX_LONG];
    int  opcion = 1;
    int  N      = 100000;

    if (argc>1) { opcion = atoi(argv[1]); } // get  first command line parameter
    if (argc>2) { N      = atoi(argv[2]); } // get second command line parameter
      
    switch(opcion) {
        case 1:
            // Secuencias de ejemplo
            strcpy(secuencia1, "actgacacgtac");
            strcpy(secuencia2, "actgacacgtac");
            
            printf("\nEjemplo 1 - Secuencias idénticas:\n");
            printf("Seq1: %s\n", secuencia1);
            printf("Seq2: %s\n", secuencia2);
            Alineamiento ej1 = alinear_secuencias_con_matriz(secuencia1, secuencia2, 1);
            mostrar_alineamiento(ej1);
            
            // Segundo ejemplo
            strcpy(secuencia1, "actgactgactgactg");
            strcpy(secuencia2, "acgtactgacgtactg");
            
            printf("\n\nEjemplo 2 - Secuencias con diferencias:\n");
            printf("Seq1: %s\n", secuencia1);
            printf("Seq2: %s\n", secuencia2);
            Alineamiento ej2 = alinear_secuencias_con_matriz(secuencia1, secuencia2, true);
            mostrar_alineamiento(ej2);
            break;
        case 2:
            printf("\n¡Bienvenido a la versión Tiled del alineador!\n");
            printf("Esta versión está en desarrollo y se lanzará próximamente.\n");
            break;                
        case 3:
            printf("\nExit\n");
            break;
            
        default:
            printf("\nNot Valid\n");
    }
    return 0;
}