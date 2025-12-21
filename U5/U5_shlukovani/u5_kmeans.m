clear; clc; close;

%% Pro generaci bodů
A = randn(23,2)*6;
B = randn(18,2)*4+[13,7];
C = randn(32,2)*8+[10,15];
D = randn(21,2)*5+[-4,9];
E = randn(23,2)*6+[-7, -17];
F = randn(18,2)*4+[23,7];
G = randn(32,2)*8+[10,-15];
H = randn(21,2)*5+[-18,9];
I = randn(32,2)*8+[32,15];
J = randn(21,2)*5+[-21,-21];
M = [A;B;C;D;E;F;G;H;I;J];
save body.mat M

load("body.mat");
figure
hold on

% plot(M(:,1), M(:,2), "g*")
%[xi,yi] = getpts;

% lineární inicializace shluků
S = [];
min_x = min(M(:,1));
min_y = min(M(:,2));
max_x = max(M(:,1));
max_y = max(M(:,2));
d_x = max_x-min_x;
d_y = max_y-min_y;

pocet_s = 10; % number of clusters

for i = 1:pocet_s   % 1 : number of clusters
    S(i,:) = [min_x + (d_x)*((i-1)/(pocet_s-1)), min_y + (d_y)*((i-1)/(pocet_s-1))];
end
V=["r." "b." "g." "m."];
V = ["r." "b." "g." "m." "c." "y." "k." "r*" "b*" "g*"];

Imax=30;
I=0;
S_prev = S - [1, 1];
while I < Imax && ~isequal(S_prev, S) % until centroids are different or iteration limit exceed
    S_prev = S;                                % store previous centroids
    L = zeros(length(M(:,1)),1);               % array of clousest centroid ids
    for i = 1:length(M)                        % for each point
        min_dist = inf;                        % sets infinite minimal dist
        closest_centroid_i = 0;                   % sets default closest centroid id
        for centroid_i = 1:length(S)              % list through centroids
            d = dist(M(i,:), S(centroid_i,:));    % calculates dist. between point and centroid
            if d < min_dist                    % if point is closer to this centroid
                min_dist = d;                  % saves new minimal dist.
                closest_centroid_i = centroid_i;     % updates closest centroid id
            end
        end
        L(i) = closest_centroid_i;                % stores closest centroid id
    end

    %%% calculates new centroids %%%
    for centroid_i = 1:length(S)                   % for each cluster
        s_x = 0;                                % sum of each x
        s_y = 0;                                % sum of each y
        count = 0;                              % count of closest points
        for i = 1:length(M)                     % for each point in cloud
            if L(i) == centroid_i                  % if this centroid is closest to the point
                s_x = s_x + M(i,1);             % add points x to sum
                s_y = s_y + M(i,2);             % add points y to sum
                count = count + 1;              % add to count of closest points
            end
        end

        if count > 0                            % if there is atleast 1 closest
            S(centroid_i, :) = [s_x/count, s_y/count];  % assigns new average x,y
        end

    end
    I = I + 1;                                  % iteration + 1
end

% vykreslení vlastní implementace
hold on
for i = 1:length(M)
    plot(M(i,1), M(i,2), V(L(i)), 'MarkerSize',12)
end
ax = plot(S(:,1), S(:,2), "kx", "MarkerSize",15, 'LineWidth',3);
title 'Implemented k-means function'
hold off

% vykreslení vestavěné metody kmeans
[idx,C] = kmeans(M,10);

figure;
axis equal
hold on
for i = 1:length(M)
    plot(M(i,1), M(i,2),V(idx(i)),'MarkerSize',12)
end
bx = plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3);
title 'Matlab k-means function'
hold off

% uložení výstupů
saveas(bx, "matlab3.png")
saveas(ax, "implemented3.png")

% funkce pro výpočet vzdálenosti
function d = dist(point1, point2)
d = sqrt((point1(1)-point2(1))^2 + (point1(2)-point2(2))^2);
end