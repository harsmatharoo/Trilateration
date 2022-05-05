% Beacon Locations
b1 = [0 0.5 0];
b2 = [3.5 0 1];
b3 = [0 1.5 0];
b4 = [4.5 0 0];
b5 = [0 2.5 3.5];
b6 = [1.5 2.5 0];
b7 = [3 0 2];
b8 = [4 4 4];
Barray = [b1; b2; b3; b4; b5; b6; b7; b8];

% convert to km/s
vcty_x = 40/1000; 
vcty_y = 20/1000; 
vcty_z = 5/1000; 

%%%%%%%%%%%%%Actual distances%%%%%%%%%%%%% - the real ones as given
Interval = [0 : 20 : 100]'; 

% x position
x1 = vcty_x.*Interval; % distance
x23 = x1(6).*(ones([10 1])); % no movement and keep previous x position 
xpos = vertcat(x1, x23); %concatenate vectors

% y position
y1 = zeros([5 1]); %no movement
y2 = vcty_y.*Interval; 
y3 = y2(6).*(ones([5 1])); 
ypos = vertcat(y1, y2, y3);

% z position
z12 = zeros([10 1]); %no movement
z3 = vcty_z.*Interval; 
zpos = vertcat(z12, z3); 

XYZ = horzcat(xpos, ypos, zpos); %Concatenate x, y and z vectors horizontally

r = []; % actual distances from beacons to target
for s = [1:1:16] 
    %copy x,y,z pos 8 times for each time interval:
    posXArr = XYZ(s,1).*ones([8 1]); 
    posYArr = XYZ(s,2).*ones([8 1]);
    posZArr = XYZ(s,3).*ones([8 1]);
    %calculate distance from all 8 beacons at each time interval
    rNew = sqrt((posXArr-Barray(1:8,1)).^2 + (posYArr-Barray(1:8,2)).^2 + (posZArr-Barray(1:8,3)).^2);
    r = vertcat(r, rNew'); % Concatenate vectors vertically
end
r 

%%%%%%%%%%%%%Estimated distances%%%%%%%%%%%%%

% Used my student number to get the seeding
rng(400185871); %random seed
stdDev = 2/1000; %standard dev in km
mean = 0;
res_array = [];
r_i = 0;
for m = [1:1:16]
    r_i = r(m, :) + stdDev.*randn(1,8) + mean;  % add noises to actual distance
    res_array = vertcat(res_array, r_i);
end

disp("Now displaying after the contribution of the added noise below:")
res_array

%%%%%%%%%%%%%Coordinates%%%%%%%%%%%%%
% temporary matrix below
Transp_mat = [0 0 0; 0 0 0; 0 0 0];
Transp_mat_New = [];
R = [0; 0; 0];
R_arr = [];

Ts = 16; % 16 points within 300 seconds

for time = [1:1:Ts]
    for k = [1:1:100] % 100 iterations
        % reset to null/zero
        Transp_mat = [0 0 0; 0 0 0; 0 0 0];
        Transp_mat_New = [];
        JF = [0;0;0];
        for i = [1:1:8]  % 8 beacons
            % fi for the ith beacon
            fi = sqrt((R(1, 1)-Barray(i,1)).^2 + (R(2, 1)-Barray(i,2)).^2 + (R(3, 1)-Barray(i,3)).^2) - res_array(time, i);
            
            % compute the J'J matrix elements seperately
            JJ_11 = (R(1, 1)-Barray(i,1))^2 / (fi + res_array(time, i))^2;
            JJ_12 = ((R(1, 1)-Barray(i,1)) * (R(2, 1)-Barray(i,2))) / (fi + res_array(time, i))^2;
            JJ_13 = ((R(1, 1)-Barray(i,1)) * (R(3, 1)-Barray(i,3))) / (fi + res_array(time, i))^2;
            JJ_21 = JJ_12;
            JJ_22 = (R(2, 1)-Barray(i,2))^2 / (fi + res_array(time, i))^2;
            JJ_23 = ((R(2, 1)-Barray(i,2)) * (R(3, 1)-Barray(i,3))) / (fi + res_array(time, i))^2;
            JJ_31 = JJ_13;
            JJ_32 = JJ_23;
            JJ_33 = (R(3, 1)-Barray(i,3))^2 / (fi + res_array(time, i))^2;

            % form the J'J matrix 
            Transp_mat_New = horzcat(JJ_11, JJ_12, JJ_13);
            Transp_mat_New = vertcat(Transp_mat_New, horzcat(JJ_21, JJ_22, JJ_23));
            Transp_mat_New = vertcat(Transp_mat_New, horzcat(JJ_31, JJ_32, JJ_33));
            Transp_mat = Transp_mat + Transp_mat_New; % sum elements 8 times
            
            % compute the J'F matrix elements seperately
            JF_11 = ((R(1, 1)-Barray(i,1))*fi) / (fi + res_array(time, i));
            JF_21 = ((R(2, 1)-Barray(i,2))*fi) / (fi + res_array(time, i));
            JF_31 = ((R(3, 1)-Barray(i,3))*fi) / (fi + res_array(time, i));

            % form the J'F matrix and sum 8 times
            JF = JF + vertcat(JF_11, JF_21, JF_31);
        end
        R = vertcat(R - (inv(Transp_mat)) * JF); % Apply the iterative formula
    end
    R_arr = vertcat(R_arr, R'); % Concatenate vectors of coordinates vertically
end
R_arr

plot3 (R_arr(6:11,1), R_arr(6:11,2),R_arr(6:11,3),'b')
hold on
plot3 (XYZ(6:11,1), XYZ(6:11,2),XYZ(6:11,3),'r')
xlabel('x axis');
ylabel('y axis');
zlabel('z axis');
legend('Estimated', 'Actual');
grid on;
