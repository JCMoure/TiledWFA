#!/usr/bin/env python3
from PIL import Image, ImageDraw, ImageFont
import os

def texto_a_imagen(archivo_txt, archivo_salida="out.png", escala=1.0):
    """
    Convierte un archivo de texto con formato de matriz a una imagen PNG
    
    Args:
        archivo_txt: Ruta al archivo de texto
        archivo_salida: Ruta de salida para la imagen
        escala: Factor de escala (1.0 = tamaño normal, mayor = más zoom)
    """
    
    # Leer el archivo de texto
    with open(archivo_txt, 'r', encoding='utf-8') as f:
        lineas = f.readlines()
    
    # Configuración de fuente
    try:
        # Intentar usar una fuente del sistema
        font_size = int(14 * escala)
        font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationMono-Regular.ttf", font_size)
    except:
        try:
            font = ImageFont.truetype("Courier", font_size)
        except:
            # Fallback a fuente por defecto
            font = ImageFont.load_default()
            font_size = 12
    
    # Calcular dimensiones del texto
    # Usar un contexto temporal para medir
    temp_img = Image.new('RGB', (1, 1))
    temp_draw = ImageDraw.Draw(temp_img)
    
    # Medir el tamaño de cada celda
    sample_text = " *999*"
    bbox = temp_draw.textbbox((0, 0), sample_text, font=font)
    cell_width = bbox[2] - bbox[0] + 10
    cell_height = bbox[3] - bbox[1] + 8
    
    # Calcular dimensiones de la imagen
    max_line_len = max(len(linea) for linea in lineas)
    
    # Contar columnas de la matriz (aprox)
    matrix_start = 0
    for i, linea in enumerate(lineas):
        if '- |' in linea or 'c |' in linea or any(c.isalpha() and '|' in linea for c in linea):
            matrix_start = i
            break
    
    # Determinar número de columnas y filas
    num_filas = 0
    num_columnas = 0
    for linea in lineas[matrix_start:]:
        if '-----' in linea or '|' not in linea:
            continue
        if '|' in linea and '=' not in linea:
            num_filas += 1
            # Contar valores en la línea
            parts = linea.split('|')
            if len(parts) > 1:
                valores = parts[1].strip().split()
                num_columnas = max(num_columnas, len(valores))
    
    # Calcular dimensiones de la imagen con padding
    padding = 20 * escala
    width = int(padding * 2 + (num_columnas + 1) * cell_width + 100)
    height = int(padding * 2 + (num_filas + 1) * cell_height + 100)
    
    # Crear imagen
    img = Image.new('RGB', (width, height), 'white')
    draw = ImageDraw.Draw(img)
    
    # Colores
    colores = {
        'normal': (0, 0, 0),
        'rojo': (255, 0, 0),
        'fondo_rojo': (255, 240, 240),
        'gris': (200, 200, 200),
        'azul_claro': (240, 248, 255),
    }
    
    # Dibujar el texto línea por línea
    y = padding
    for idx, linea in enumerate(lineas):
        # Para la matriz, procesar cada línea
        if '-----' in linea:
            # Línea separadora
            draw.line([(padding, y + cell_height/2), 
                       (width - padding, y + cell_height/2)], 
                      fill=colores['gris'], width=2)
            y += cell_height/2
            continue
        
        # Dibujar línea de texto
        x = padding
        i = 0
        while i < len(linea):
            # Detectar celdas especiales con *
            if i < len(linea) - 5 and linea[i:i+2] == ' *':
                # Celda especial (en rojo)
                end = i + 2
                while end < len(linea) and linea[end] != '*':
                    end += 1
                if end < len(linea):
                    texto_celda = linea[i+2:end+1]  # Incluye el * final
                    draw.rectangle([x-2, y-2, x+cell_width+2, y+cell_height+2], 
                                 fill=colores['azul_claro'])
                    draw.text((x, y), texto_celda, fill=colores['rojo'], font=font)
                    x += cell_width
                    i = end + 1
                    continue
            
            # Texto normal
            if linea[i] != ' ' or (i > 0 and linea[i-1] == ' '):
                # Encontrar el final de la palabra
                end = i
                while end < len(linea) and linea[end] != ' ':
                    end += 1
                if end > i:
                    palabra = linea[i:end]
                    draw.text((x, y), palabra, fill=colores['normal'], font=font)
                    # Calcular ancho aproximado
                    bbox = draw.textbbox((x, y), palabra, font=font)
                    x += bbox[2] - bbox[0] + 5
                    i = end
                    continue
            i += 1
        
        y += cell_height
    
    # Guardar imagen
    img.save(archivo_salida)
    print(f"Imagen guardada como {archivo_salida}")
    print(f"Dimensiones: {width}x{height} píxeles")

if __name__ == "__main__":
    # Verificar si existe el archivo
    if not os.path.exists("matriz.txt"):
        print("Error: No se encuentra el archivo matriz.txt")
        print("Asegúrate de ejecutar primero el programa C que genera el archivo")
        exit(1)
    
    # Convertir a PNG
    texto_a_imagen("matriz.txt", "out.png", escala=1.0)
    
    # Opcional: También generar una versión con más zoom
    texto_a_imagen("matriz.txt", "out_zoom.png", escala=1.5)
    print("\nVersión con zoom generada: out_zoom.png")