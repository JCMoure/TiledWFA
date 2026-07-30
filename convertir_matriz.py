#!/usr/bin/env python3
from PIL import Image, ImageDraw, ImageFont
import os
import sys

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
    font_size = int(14 * escala)
    try:
        # Intentar usar una fuente del sistema
        font = ImageFont.truetype("/usr/share/fonts/truetype/liberation/LiberationMono-Regular.ttf", font_size)
    except:
        try:
            font = ImageFont.truetype("Courier", font_size)
        except:
            # Fallback a fuente por defecto
            font = ImageFont.load_default()
            font_size = 12
    
    # Función para medir texto (compatible con versiones antiguas de Pillow)
    def get_text_size(text, font):
        """Obtiene el tamaño de un texto usando métodos compatibles"""
        temp_img = Image.new('RGB', (1, 1))
        temp_draw = ImageDraw.Draw(temp_img)
        
        # Método 1: Usar textbbox (nuevo)
        try:
            bbox = temp_draw.textbbox((0, 0), text, font=font)
            return bbox[2] - bbox[0], bbox[3] - bbox[1]
        except AttributeError:
            # Método 2: Usar textsize (antiguo, pero compatible)
            try:
                size = temp_draw.textsize(text, font=font)
                return size[0], size[1]
            except AttributeError:
                # Método 3: Estimar tamaño
                return len(text) * font_size * 0.6, font_size
    
    # Calcular el tamaño de cada celda
    sample_text = " *999*"
    text_width, text_height = get_text_size(sample_text, font)
    cell_width = int(text_width + 10 * escala)
    cell_height = int(text_height + 8 * escala)
    
    # Determinar estructura de la matriz
    # Encontrar dónde comienza la matriz (línea con formato de tabla)
    matrix_start = 0
    for i, linea in enumerate(lineas):
        if '- |' in linea or ('|' in linea and i > 2):
            matrix_start = i
            break
    
    # Contar filas y columnas de la matriz
    num_filas = 0
    num_columnas = 0
    for linea in lineas[matrix_start:]:
        if '-----' in linea or '|' not in linea:
            continue
        if '|' in linea and '=' not in linea and '---' not in linea:
            num_filas += 1
            # Contar valores en la línea (separados por |)
            partes = linea.split('|')
            if len(partes) > 1:
                # Contar elementos en la parte derecha
                valores = partes[1].strip().split()
                num_columnas = max(num_columnas, len(valores))
    
    # Si no se detectaron columnas, usar valor predeterminado
    if num_columnas == 0:
        num_columnas = 10  # Valor predeterminado
    
    # Calcular dimensiones de la imagen con padding
    padding = int(20 * escala)
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
        'gris': (180, 180, 180),
        'azul_claro': (230, 240, 255),
        'separador': (200, 200, 200),
    }
    
    # Función para dibujar texto con compatibilidad
    def draw_text(draw, pos, text, fill, font):
        try:
            draw.text(pos, text, fill=fill, font=font)
        except:
            # Fallback si no se puede usar la fuente
            draw.text(pos, text, fill=fill)
    
    # Dibujar el texto línea por línea
    y = padding
    fila_actual = 0
    
    for idx, linea in enumerate(lineas):
        linea = linea.rstrip('\n')
        
        # Saltar líneas vacías
        if not linea.strip():
            y += cell_height / 2
            continue
        
        # Línea separadora de la matriz (-----)
        if '-----' in linea:
            y_sep = y + cell_height/2
            draw.line([(padding, y_sep), (width - padding, y_sep)], 
                      fill=colores['separador'], width=2)
            y += cell_height/2
            continue
        
        # Procesar líneas de la matriz (contienen '|')
        if '|' in linea and idx >= matrix_start:
            # Dividir la línea en partes
            partes = linea.split('|')
            
            # Parte izquierda (encabezado de fila)
            if len(partes) > 0 and partes[0].strip():
                x = padding + 5
                texto_header = partes[0].strip()
                draw_text(draw, (x, y + 5), texto_header, colores['normal'], font)
            
            # Parte derecha (valores de la matriz)
            if len(partes) > 1:
                # Extraer valores, incluyendo los marcadores *
                valores = partes[1].strip().split()
                x = padding + cell_width + 10
                
                for valor in valores:
                    # Verificar si es un valor especial (con *)
                    es_especial = '*' in valor
                    
                    # Limpiar el valor para mostrar
                    valor_limpio = valor.replace('*', '')
                    
                    # Dibujar fondo para celdas especiales
                    if es_especial:
                        draw.rectangle([x-3, y-2, x+cell_width-5, y+cell_height-2], 
                                     fill=colores['azul_claro'])
                    
                    # Dibujar el valor
                    color = colores['rojo'] if es_especial else colores['normal']
                    draw_text(draw, (x + 5, y + 5), valor_limpio, color, font)
                    
                    # Dibujar bordes de celda (excepto para la última)
                    draw.rectangle([x, y, x+cell_width-2, y+cell_height-2], 
                                 outline=colores['gris'], width=1)
                    
                    x += cell_width
            
            y += cell_height
            fila_actual += 1
        else:
            # Texto normal (títulos, etc.)
            # Para encabezados, centrar
            if idx < matrix_start:
                draw_text(draw, (padding, y + 5), linea, colores['normal'], font)
                y += cell_height
            else:
                # Otro texto fuera de la matriz
                if linea.strip():
                    draw_text(draw, (padding, y + 5), linea, colores['normal'], font)
                    y += cell_height / 2
    
    # Guardar imagen
    img.save(archivo_salida)
    print(f"✅ Imagen guardada como {archivo_salida}")
    print(f"📐 Dimensiones: {width}x{height} píxeles")
    print(f"📊 Filas: {num_filas}, Columnas: {num_columnas}")

if __name__ == "__main__":
    # Verificar si existe el archivo
    if not os.path.exists("matriz.txt"):
        print("❌ Error: No se encuentra el archivo matriz.txt")
        print("Asegúrate de ejecutar primero el programa C que genera el archivo")
        sys.exit(1)
    
    # Verificar versión de Pillow
    try:
        from PIL import __version__ as pillow_version
        print(f"📦 Pillow versión: {pillow_version}")
    except:
        pass
    
    # Convertir a PNG
    print("🔄 Convirtiendo matriz.txt a imagen...")
    texto_a_imagen("matriz.txt", "out.png", escala=1.0)
    
    # Generar versión con zoom
    print("\n🔍 Generando versión con zoom...")
    texto_a_imagen("matriz.txt", "out_zoom.png", escala=1.8)
    
    print("\n✨ Archivos generados:")
    print("   - out.png (tamaño normal)")
    print("   - out_zoom.png (con zoom)")
    
    # Intentar abrir la imagen
    try:
        if sys.platform == 'darwin':  # macOS
            os.system(f"open out.png")
        elif sys.platform == 'linux':
            # Intentar diferentes visores de imágenes
            for viewer in ['display', 'eog', 'xdg-open', 'gimp']:
                if os.system(f"which {viewer} > /dev/null 2>&1") == 0:
                    os.system(f"{viewer} out.png &")
                    break
        elif sys.platform == 'win32':  # Windows
            os.system(f"start out.png")
    except:
        print("💡 Puedes abrir out.png manualmente con tu visor de imágenes favorito")