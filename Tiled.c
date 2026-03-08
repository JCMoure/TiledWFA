#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX_LONG 100
#define MATCH 2
#define MISMATCH -1
#define GAP -2

// Estructura para almacenar el resultado del alineamiento
typedef struct {
    char seq1_alineada[MAX_LONG * 2];
    char seq2_alineada[MAX_LONG * 2];
    int puntuacion;
} Alineamiento;

// Función para calcular el máximo de tres números
int maximo(int a, int b, int c) {
    int max = a;
    if (b > max) max = b;
    if (c > max) max = c;
    return max;
}

// Función para obtener la puntuación de match/mismatch
int obtener_puntuacion(char a, char b) {
    // Convertir a minúsculas para comparación
    a = tolower(a);
    b = tolower(b);
    
    if (a == b) {
        return MATCH;
    } else {
        return MISMATCH;
    }
}

// Función principal de alineamiento (Needleman-Wunsch)
Alineamiento alinear_secuencias(char *seq1, char *seq2) {
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
            
            matriz[i][j] = maximo(match, del, ins);
            
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
    
    // Trazado para obtener el alineamiento
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
    
    // Invertir las secuencias (las construimos al revés)
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

// Función para validar secuencias de ADN
int validar_secuencia(char *secuencia) {
    for (int i = 0; i < strlen(secuencia); i++) {
        char c = tolower(secuencia[i]);
        if (c != 'a' && c != 'c' && c != 'g' && c != 't') {
            return 0; // Carácter no válido
        }
    }
    return 1; // Secuencia válida
}

// Función para limpiar el buffer de entrada
void limpiar_buffer() {
    int c;
    while ((c = getchar()) != '\n' && c != EOF);
}

int main() {
    char secuencia1[MAX_LONG];
    char secuencia2[MAX_LONG];
    int opcion;
      
    do {
        printf("\n--- MENÚ PRINCIPAL ---\n");
        printf("1. Usar secuencias de ejemplo\n");
        printf("2. Tiled Version\n");
        printf("3. Salir\n");
        printf("Selecciona una opción: ");
        
        scanf("%d", &opcion);
        limpiar_buffer();
        
        switch(opcion) {
            case 1:
                // Secuencias de ejemplo
                strcpy(secuencia1, "actgac");
                strcpy(secuencia2, "actgac");
                
                printf("\nEjemplo 1 - Secuencias idénticas:\n");
                printf("Seq1: %s\n", secuencia1);
                printf("Seq2: %s\n", secuencia2);
                Alineamiento ej1 = alinear_secuencias(secuencia1, secuencia2);
                mostrar_alineamiento(ej1);
                
                // Segundo ejemplo
                strcpy(secuencia1, "actg");
                strcpy(secuencia2, "acgt");
                
                printf("\n\nEjemplo 2 - Secuencias con diferencias:\n");
                printf("Seq1: %s\n", secuencia1);
                printf("Seq2: %s\n", secuencia2);
                Alineamiento ej2 = alinear_secuencias(secuencia1, secuencia2);
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
        
    } while (opcion != 3);
    
    return 0;
}