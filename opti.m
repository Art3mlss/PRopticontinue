%% S7 - Projet OPTI - RICALENS BARON

%% Question 1)

% Paramètres
R = 1.5;  % Rayon du cercle
n = 30;  % Nombre de points

% Génération des points (x_i, y_i) autour du cercle de centre (0,0)
theta = linspace(0, 2*pi, n);

% Définition des intervalles de cx et cy
cx_vals1 = linspace(-1, 1, 100);
cy_vals1 = linspace(-1, 2, 100);
cx_vals2 = linspace(-1, 4, 100);
cy_vals2 = linspace(-1, 4, 100);

% Calcul de la fonction de coût pour les deux intervalles
[Cx1, Cy1] = meshgrid(cx_vals1, cy_vals1);
[Cx2, Cy2] = meshgrid(cx_vals2, cy_vals2);

% Fonction de coût CTLS pour les deux grilles
CTLS1 = zeros(size(Cx1));
CTLS2 = zeros(size(Cx2));

for i = 1:n
    % Pour chaque point (x_i, y_i), calcul de l'écart entre Di et R
    CTLS1 = CTLS1 + (sqrt((xi(i) - Cx1).^2 + (yi(i) - Cy1).^2) - R).^2;
    CTLS2 = CTLS2 + (sqrt((xi(i) - Cx2).^2 + (yi(i) - Cy2).^2) - R).^2;
end

% Paramètres
R = 1.5;  % Rayon du cercle
n = 30;  % Nombre de points

% Génération des points (x_i, y_i) autour du cercle de centre (0,0)
theta = linspace(0, 2*pi, n);
x = R * cos(theta) + 0.1*randn(1, n); % ajout de bruit pour simuler les erreurs de mesure
y = R * sin(theta) + 0.1*randn(1, n);

% Définition des intervalles de cx et cy
cx_vals1 = linspace(-1, 1, 100);
cy_vals1 = linspace(-1, 2, 100);
cx_vals2 = linspace(-1, 4, 100);
cy_vals2 = linspace(-1, 4, 100);

% Calcul de la fonction de coût pour les deux intervalles
[Cx1, Cy1] = meshgrid(cx_vals1, cy_vals1);
[Cx2, Cy2] = meshgrid(cx_vals2, cy_vals2);

% Fonction de coût CTLS pour les deux grilles
CTLS1 = zeros(size(Cx1));
CTLS2 = zeros(size(Cx2));

for i = 1:n
    % Pour chaque point (x_i, y_i), calcul de l'écart entre Di et R
    CTLS1 = CTLS1 + (sqrt((xi(i) - Cx1).^2 + (yi(i) - Cy1).^2) - R).^2;
    CTLS2 = CTLS2 + (sqrt((xi(i) - Cx2).^2 + (yi(i) - Cy2).^2) - R).^2;
end

% Représentation graphique
figure;

% 1ère sous-figure: Vue 3D sur l'intervalle [-1, 1] × [-1, 2]
subplot(2,2,1);
surf(Cx1, Cy1, CTLS1);  % Affichage 3D
title('CTLS en 3D sur [-1, 1] × [-1, 2]');
xlabel('cx'); ylabel('cy'); zlabel('CTLS');
shading interp;  % Pour lisser les couleurs
view(3);  % Vue 3D

% 2ème sous-figure: Vue en 2D (contour) sur l'intervalle [-1, 1] × [-1, 2]
subplot(2,2,2);
contourf(Cx1, Cy1, CTLS1, 20);  % Affichage 2D en niveaux de contour
title('CTLS en 2D (contour) sur [-1, 1] × [-1, 2]');
xlabel('cx'); ylabel('cy'); zlabel('CTLS');
colorbar;

% 3ème sous-figure: Vue 3D sur l'intervalle [-1, 4] × [-1, 4]
subplot(2,2,3);
surf(Cx2, Cy2, CTLS2);  % Affichage 3D
title('CTLS en 3D sur [-1, 4] × [-1, 4]');
xlabel('cx'); ylabel('cy'); zlabel('CTLS');
shading interp;  % Pour lisser les couleurs
view(3);  % Vue 3D

% 4ème sous-figure: Vue en 2D (contour) sur l'intervalle [-1, 4] × [-1, 4]
subplot(2,2,4);
contourf(Cx2, Cy2, CTLS2, 20);  % Affichage 2D en niveaux de contour
title('CTLS en 2D (contour) sur [-1, 4] × [-1, 4]');
xlabel('cx'); ylabel('cy'); zlabel('CTLS');
colorbar;

%% Question 2)

precision = 1e-2;  % Précision désirée

% Grilles de recherche pour les deux domaines
cx_vals1 = -1:precision:1;
cy_vals1 = -1:precision:2;
cx_vals2 = -1:precision:4;
cy_vals2 = -1:precision:4;

% Initialisation des matrices pour stocker les coûts
CTLS1 = zeros(length(cy_vals1), length(cx_vals1));
CTLS2 = zeros(length(cy_vals2), length(cx_vals2));

% Fonction de coût sur le premier domaine [-1, 1] x [-1, 2]
for i = 1:n
    for cx_idx = 1:length(cx_vals1)
        for cy_idx = 1:length(cy_vals1)
            cx = cx_vals1(cx_idx);
            cy = cy_vals1(cy_idx);
            % Calcul de la fonction de coût
            CTLS1(cy_idx, cx_idx) = CTLS1(cy_idx, cx_idx) + (sqrt((xi(i) - cx)^2 + (yi(i) - cy)^2) - R)^2;
        end
    end
end

% Fonction de coût sur le second domaine [-1, 4] x [-1, 4]
for i = 1:n
    for cx_idx = 1:length(cx_vals2)
        for cy_idx = 1:length(cy_vals2)
            cx = cx_vals2(cx_idx);
            cy = cy_vals2(cy_idx);
            % Calcul de la fonction de coût
            CTLS2(cy_idx, cx_idx) = CTLS2(cy_idx, cx_idx) + (sqrt((xi(i) - cx)^2 + (yi(i) - cy)^2) - R)^2;
        end
    end
end

% Recherche des minima pour les deux domaines
[min_CTLS1, idx1] = min(CTLS1(:));
[min_CTLS2, idx2] = min(CTLS2(:));

% Conversion des indices en coordonnées cx, cy
[cy_min_idx1, cx_min_idx1] = ind2sub(size(CTLS1), idx1);
[cy_min_idx2, cx_min_idx2] = ind2sub(size(CTLS2), idx2);

cx_min1 = cx_vals1(cx_min_idx1);
cy_min1 = cy_vals1(cy_min_idx1);

cx_min2 = cx_vals2(cx_min_idx2);
cy_min2 = cy_vals2(cy_min_idx2);

% Affichage des résultats
disp(['Minimium dans le premier domaine: cx = ', num2str(cx_min1), ', cy = ', num2str(cy_min1)]);
disp(['Minimium dans le second domaine: cx = ', num2str(cx_min2), ', cy = ', num2str(cy_min2)]);

% Représentation des points et des cercles trouvés
figure;
hold on;
plot(xi, yi, 'b*');  % Nuage de points
viscircles([cx_min1, cy_min1], R, 'Color', 'r');  % Cercle trouvé pour le 1er domaine
viscircles([cx_min2, cy_min2], R, 'Color', 'g');  % Cercle trouvé pour le 2e domaine
axis equal;
hold off;

% Cela montre que notre fonction cout est merdique, faut la changer en fait
% pour avoir