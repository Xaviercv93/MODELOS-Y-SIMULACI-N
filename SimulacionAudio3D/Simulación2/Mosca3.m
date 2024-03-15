% Limpiar el espacio de trabajo
close all
clear
clc

%Cargar el archivo de audio
[audio, Fs] = audioread('moscaShort.wav');

%numero de muestras en el vector audio
longitud_audio= length(audio);
%frecuencia de muestreo
Fs1 = 44100;

% Leer información del archivo de audio
info = audioinfo('moscaShort.wav');
% Obtener el número de canales
num_canales = info.NumChannels;
% Mostrar el número de canales
disp(['Número de canales: ', num2str(num_canales),newline]);

%duracion original del audio
longitud_audio=longitud_audio/Fs1;
disp(['Muestras del audio modificado: ', num2str(length(audio))])
disp(['Longitud del audio: ', num2str(longitud_audio)]);
%----------------DURACION AUDIO------------
% Limitar la duración total del audio a 5 segundos
duracion_total_deseada = 5; % Duración deseada en segundos
num_muestras_deseadas = duracion_total_deseada * Fs1; % Convertir duración a número de muestras
audio = audio(1:min(length(audio), num_muestras_deseadas), :); % Tomar solo las primeras num_muestras_deseadas muestras
disp(['Muestras del audio modificado: ', num2str(length(audio))])
disp(['Duracion del audio modificado: ', num2str(length(audio)/Fs1)]);


%Duración de cada bloque en segundos
duracion_segmento = 0.5; % Por ejemplo, 0.1 segundos

%Número de muestras por bloque
muestras_por_segmentos = duracion_segmento * Fs1;

%Calcular el número total de bloques
num_segmentos = ceil(length(audio) / muestras_por_segmentos);
fprintf('\nduracion_bloque: %.1f\nmuestras_por_bloque:    %d\nnum_bloques:  %d\n', duracion_segmento, muestras_por_segmentos,num_segmentos);

%Calcular la ventana de crossfade exponencial decreciente y creciente
crossfade_length = 512; % Longitud del crossfade en muestras
crossfade_samples = min(crossfade_length, muestras_por_segmentos); % Asegurarse de que el crossfade no exceda la longitud del bloque

% Funcion exponencial creciente y decreciente
fade_exp_out = exp(-(0:crossfade_samples-1)/(crossfade_samples-1) * 1.5); 
fade_exp_in = exp((0:crossfade_samples-1)/(crossfade_samples-1) * 1.5); 

%Filtro de la Resuesta Impulsiva para la convolución
H = fir1(100, 0.1); % Ejemplo de filtro de la Respuesta impulsiva

%Dividir el audio en canales (izquierdo y derecho)
audio_L = audio(:, 1);
audio_R = audio(:, 2);

%Parámetros de simulación de azimut
azimut_max = 360; % Azimut máximo

for segmentos = 1:num_segmentos
    %Calcular el azimut para este segmento, mas segmentos aplica mas azimut
    azimut = (segmentos - 1) * (azimut_max / num_segmentos); 
    
    %Debido a la posicion de los oidos se calcular la ganancia para cada canal (izquierdo y derecho) basada en el azimut
    ganancia_izquierdo = cosd(azimut+90);
    ganancia_derecho = cosd(azimut);
    
    %En el segmento actual calcular el índice de inicio y final
    inicio_index = (segmentos - 1) * muestras_por_segmentos + 1;
    fin_index = min(segmentos * muestras_por_segmentos, length(audio));
    
    %Extraer segmento del cada canal
    segmento_actual_L = audio_L(inicio_index:fin_index);
    segmento_actual_R = audio_R(inicio_index:fin_index);
    
    %Aplicar crossfade exponencial decreciente a segmento de canal izquierdo
    if segmentos > 1
        fade_out_section = fade_exp_out(end-crossfade_samples+1:end);
        segmento_actual_L(1:crossfade_samples) = segmento_actual_L(1:crossfade_samples) .* fade_out_section';
    end
    
    %Aplicar crossfade exponencial creciente a segmento de canal derecho
    if segmentos < num_segmentos
        fade_in_section = fade_exp_in(end - crossfade_samples + 1 : end);
        segmento_actual_R(end - crossfade_samples + 1 : end) = segmento_actual_R(end - crossfade_samples + 1 : end) .* fade_in_section';
    end
    
    %Añadir ganancia dado el azimut en cada canal
    segmento_actual_L = segmento_actual_L * ganancia_izquierdo;
    segmento_actual_R = segmento_actual_R * ganancia_derecho;
    
    %Convolucionar de filtro FIR y segmento
    segmento_actual_L = conv(segmento_actual_L, H, 'same');
    segmento_actual_R = conv(segmento_actual_R, H, 'same');
    
    %Agregar segmento actual a salida
    audio_L(inicio_index:fin_index) = segmento_actual_L;
    audio_R(inicio_index:fin_index) = segmento_actual_R;
end

%Combinar canales
output_audio = [audio_L, audio_R];

sound(output_audio, Fs1);

pause(length(output_audio) / Fs1);

%Graficar
t = (0:length(output_audio)-1) / Fs1;
figure;
subplot(2, 1, 1);
plot(t, audio_L, 'b');
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Canal Izquierdo');
subplot(2, 1, 2);
plot(t, audio_R, 'r');
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Canal Derecho');