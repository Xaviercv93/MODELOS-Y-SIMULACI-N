

[audio, Fs] = audioread('chord.wav');

% duración de los segmentos en segundos
duracion_segmento = 10;

% Cantidad de muestras por segmento
muestras_por_segmento = duracion_segmento * Fs;

% Cantidad total de segmentos
cantidad_segmentos = ceil(length(audio) / muestras_por_segmento);

%arreglo de segmentos
segmentos_audio = zeros(muestras_por_segmento, cantidad_segmentos);

% Dividir el audio en segmentos
for i = 1:cantidad_segmentos
    inicio = (i - 1) * muestras_por_segmento + 1;
    fin = min(i * muestras_por_segmento, length(audio));
    segmento_actual = audio(inicio:fin, :);
    segmentos_audio(1:length(segmento_actual), i) = segmento_actual(i);
end

%H = interphrir(90, 0);

Ts = 1 / Fs;
ele = 0;
crossfade_duracion = 0.1;  % Duración del crossfade en segundos

for a = 0:90:360
    for i = 1:size(segmentos_audio, 2)-1
        % Interpolación
        H=interphrir(a,ele);
        % Interpolación lineal para el crossfade
        t_crossfade = linspace(0, 1, crossfade_duracion * Fs);
        crossfade = 1 - t_crossfade;

        % Se aplica el crossfade a los segmentos
        sd(:, 1) = conv(segmentos_audio(:, i), H(:, 1)) * crossfade' + ...
                    conv(segmentos_audio(:, i+1), H(:, 1)) * t_crossfade';
        
        sd(:, 2) = conv(segmentos_audio(:, i), H(:, 2)) * crossfade' + ...
                    conv(segmentos_audio(:, i+1), H(:, 2)) * t_crossfade';
        
        % se normaliza la señal resultante
        sd = normalize(sd, 'range');
        
        % Reproducir la señal
        sound(sd, Fs);
        pause(0.3)


    end
end
