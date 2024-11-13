%% S7 - Projet OPTI - RICALENS BARON

%% Question 1)

% Paramètres
R = 1.5;  % Rayon du cercle
n = 30;  % Nombre de points

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

%% Question 4)

% Fonction de calcul
dbtype("gradient_CTLS.m");

% Paramètres pour la vérification
cx_test = 0.5;
cy_test = 1.0;
R = 1.5;
epsilon = 1e-5;  % Petite perturbation pour les différences finies

% Points de mesure
theta = linspace(0, 2*pi, 100);
x = R * cos(theta) + 0.05 * randn(1, 100);  % Bruit ajouté pour simuler les erreurs
y = R * sin(theta) + 0.05 * randn(1, 100);

function cost = CTLS(cx, cy)
    global xi yi R;
    cost = 0;
    
    for i = 1:size(xi)
        Di = sqrt((xi(i) - cx)^2 + (yi(i) - cy)^2);
        cost = cost + (Di - R)^2;
    end
end

% Gradient analytique
grad_analytique = gradient_CTLS(cx_test, cy_test);

% Gradient numérique par différences finies
cost_plus_cx = CTLS(cx_test + epsilon, cy_test);
cost_minus_cx = CTLS(cx_test - epsilon, cy_test);
grad_num_cx = (cost_plus_cx - cost_minus_cx) / (2 * epsilon);

cost_plus_cy = CTLS(cx_test, cy_test + epsilon);
cost_minus_cy = CTLS(cx_test, cy_test - epsilon);
grad_num_cy = (cost_plus_cy - cost_minus_cy) / (2 * epsilon);

grad_numerique = [grad_num_cx; grad_num_cy];

% Affichage des résultats
disp('Gradient analytique :');
disp(grad_analytique);
disp('Gradient numérique (différences finies) :');
disp(grad_numerique);
disp('Différence entre les deux gradients :');
disp(norm(grad_analytique - grad_numerique));

%% Question 5)

% Tracer les lignes de niveaux de la fonction de coût
figure;

contourf(Cx2, Cy2, CTLS2, 20);  % Affichage 2D en niveaux de contour
title('CTLS en 2D (contour) sur [-1, 4] × [-1, 4]');
xlabel('cx'); ylabel('cy'); zlabel('CTLS');
colorbar;
hold on;

% Calcul de la fonction de coût et du gradient sur la grille
for i = 1:size(Cx2)
    for j=1:size(Cy2)
        cx = Cx2(i, j);
        cy = Cy2(i, j);
        
        % Calcul du gradient
        grad = gradient_CTLS(cx, cy);
    end
end

% Représenter le champ de vecteurs du gradient
quiver(CX, CY, grad_x, grad_y, 'r');  % Vecteurs en rouge

% Ajustements de l'affichage
axis equal;  % Repère orthonormé
title('Champ de vecteurs du gradient et lignes de niveau de CTLS');
xlabel('cx');
ylabel('cy');
hold off;
