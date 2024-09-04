clc;
clear all;

even=input("no of lamina FOR EVEN (1) FOR ODD (0): ");
n_lamina = input('Enter the number of lamina: ');
disp('==========================================================')

%UNCOMMENT THE PART IF YOU WANT TO ASK THE INPUT FROM USER

%======================================================================================================
%Loop to input each angle
% for i = 1:n_lamina
%     angle_degrees = input(['Enter angle ' num2str(i) ' for '  ' lamina ' num2str(i) ': ']);
%      % Store the angle in the array
%     theta(i) = angle_degrees;
% % end

% E1 = input(' Longitudinal stiffness of each lamina:');
% E2 =input(' Transverse stiffness of each lamina:');
% G12 =input(' Shear stiffness of each lamina:');
% sigma1_t =input(' Longitudinal tensile strength of each lamina:');
% sigma1_c = input(' Longitudinal compressive strength of each lamina:');
% sigma2_t = input(' Transverse tensile strength of each lamina:');
% sigma2_c = input('  Transverse compressive strength of each lamina:');
% tau12 = input(' Shear strength of each lamina:');
% Nx = input('  Applied load in x-direction:');
% Ny = input(' Applied load in y-direction:');
% Nxy = input(' Applied shear load:');
% Mx = input('  Applied moment in x-direction:');
% My = input('  Applied moment in y-direction:');
% Mxy = input(' Applied twisting moment:');
% dT =input(' Temperature change:');
% dC = input('  Moisture change:');
% v12 = input('poissons ratio v12 : ');
% tk=input('thickness of each lamina: ');
%====================================================================================================

%QUESTION 1  (4 LAMINA) 

theta=[0,90,90,0];
E1 = 140; % Longitudinal stiffness of each lamina
E2 = 10;  % Transverse stiffness of each lamina
G12 = 7 ; % Shear stiffness of each lamina
sigma1_t = 1400 ;  % Longitudinal tensile strength of each lamina
sigma1_c = 1400;  % Longitudinal compressive strength of each lamina
sigma2_t = 50 ;  % Transverse tensile strength of each lamina
sigma2_c = 200 ;  % Transverse compressive strength of each lamina
tau12 = 70 ;  % Shear strength of each lamina
Nx =100;  % Applied load in x-direction
Ny = 0;  % Applied load in y-direction
Nxy = 0;  % Applied shear load
Mx = 0;  % Applied moment in x-direction
My = 0;  % Applied moment in y-direction
Mxy = 0;  % Applied twisting moment
dT = 0;  % Temperature change
dC = 0;  % Moisture change
v12=0.3;
tk=0.25;%thickness of each lamina
%======================================================================================================
% 
% %QUESTION 2 (8 LAMINA)
% 
% theta=[0,45,-45,90,90,-45,45,0];
% E1 = 38.6; % Longitudinal stiffness of each lamina
% E2 = 8.27;  % Transverse stiffness of each lamina
% G12 = 4.140 ; % Shear stiffness of each lamina
% sigma1_t = 1062 ;  % Longitudinal tensile strength of each lamina
% sigma1_c = 610;  % Longitudinal compressive strength of each lamina
% sigma2_t = 31 ;  % Transverse tensile strength of each lamina
% sigma2_c = 118 ;  % Transverse compressive strength of each lamina
% tau12 = 72 ;  % Shear strength of each lamina
% Nx =100;  % Applied load in x-direction
% Ny = 0;  % Applied load in y-direction
% Nxy = 0;  % Applied shear load
% Mx = 0;  % Applied moment in x-direction
% My = 0;  % Applied moment in y-direction
% Mxy = 0;  % Applied twisting moment
% dT = 50;  % Temperature change
% dC = 0;  % Moisture change
% v12=0.28;
% tk=0.125;
% %=================================================================================================================

alpha_1=8.6*1.0e-6;% m/m/oC
alpha_2=22.1*1.0e-6;%m/m/oC
alpha1=[alpha_1;alpha_2;0];
beta_1=1.0e-4;
beta_2=1.0e-3;
beta1=[beta_1;beta_2;0];
syms Longitudinal_direction;
syms Transverse_direction
syms shear_failure;
failure_criteria=[Longitudinal_direction,Transverse_direction,shear_failure];

% Initialize an empty array to store the angles
angles = zeros(1, n_lamina);
% Total thickness
total_thickness = n_lamina*tk;
% Initialize the array to store the thicknesses
z = zeros(1, n_lamina+1);

% Loop to calculate and store the thickness of each part
for i = 1:n_lamina+1
    % Calculate the thickness for this part
    thickness = -total_thickness/2 + (i - 1) * tk;

    % Store the thickness in the array
    z(i)=thickness;
end

v21=(E2/E1)*v12;
Q11 = E1 / (1 - v12*v21);
Q12 = v12 * E2 / (1 - v12*v21);
Q21= v12 * E2 / (1 - v12*v21);
Q16 = 0; % Assuming plane stress condition
Q22 = E2/ (1 - v12*v21);
Q26 = 0; % Assuming plane stress condition
Q61=0;
Q62=0;
Q66 = G12;
Q=[Q11 Q12 Q16;Q21 Q22 Q26;Q61 Q62 Q66];

N=[Nx;Ny;Nxy;Mx;My;Mxy];

failure_LOAD=zeros(6,n_lamina);

    plyfail=[];
    unique_values = unique(theta);
    n=length(unique_values);
    for k_outer = 1:n

        A=zeros(3,3);B=zeros(3,3);D=zeros(3,3);NT=zeros(3,1);NH=zeros(3,1);MT=zeros(3,1);MH=zeros(3,1);
        for k=1:n_lamina
            if any(plyfail==theta(k))

                continue
            else

                Qbar=zeros(3,3);
                T=  [cosd(theta(k))^2, sind(theta(k))^2, 2*sind(theta(k))*cosd(theta(k));
                    sind(theta(k))^2, cosd(theta(k))^2, -2*sind(theta(k))*cosd(theta(k));
                    -sind(theta(k))*cosd(theta(k)), sind(theta(k))*cosd(theta(k)), cosd(theta(k))^2 - sind(theta(k))^2];

                alpha_x= inv(T)*alpha1;
                alphax=[alpha_x(1,1);alpha_x(2,1);2*alpha_x(3,1)];
                beta_x= inv(T)*beta1;
                betax=[beta_x(1,1);beta_x(2,1);2*beta_x(3,1)];

                Qbar_11= Q11*(cosd(theta(k)))^4 + 2*(Q12+2*Q66)*(cosd(theta(k))^2) * (sind(theta(k)))^2 + Q22*(sind(theta(k)))^4;
                Qbar_12= (Q11+Q22-4*Q66)*(cosd(theta(k)))^2 * (sind(theta(k)))^2 + Q12*((cosd(theta(k)))^4 + (sind(theta(k)))^4);
                Qbar_16= (Q11-Q12-2*Q66)*cosd(theta(k))^3 * sind(theta(k)) - (Q22-Q12-2*Q66)*cosd(theta(k))*sind(theta(k))^3;
                Qbar_21= (Q11+Q22-4*Q66)*cosd(theta(k))^2 * sind(theta(k))^2 + Q12*(cosd(theta(k))^4 + sind(theta(k))^4);
                Qbar_22= Q11*sind(theta(k))^4 + 2*(Q12+2*Q66)*cosd(theta(k))^2 * sind(theta(k))^2 + Q22*cosd(theta(k))^4;
                Qbar_26= (Q11-Q12-2*Q66)*sind(theta(k))^3 * cosd(theta(k))-(Q22-Q12-2*Q66)*cosd(theta(k))^3 * sind(theta(k));
                Qbar_61= (Q11-Q12-2*Q66)*cosd(theta(k))^3 * sind(theta(k)) - (Q22-Q12-2*Q66)*cosd(theta(k))*sind(theta(k))^3;
                Qbar_62= (Q11-Q12-2*Q66)*sind(theta(k))^3*cosd(theta(k))-(Q22-Q12-2*Q66)*cosd(theta(k))^3*sind(theta(k));
                Qbar_66= (Q11+Q22-2*Q12-2*Q66)*cosd(theta(k))^2 * sind(theta(k))^2 + Q66*(cosd(theta(k))^4+sind(theta(k))^4);

                Qbar =[Qbar_11 Qbar_12 Qbar_16; Qbar_21 Qbar_22 Qbar_26 ; Qbar_61 Qbar_62 Qbar_66];
                % Qbar= Qbar + inv(T)*Q*R*T

                A= A + Qbar.*(z(k+1)-z(k));
                B= B + 0.5.* Qbar.*(z(k+1)^2-z(k)^2);
                D= D + (1/3).* Qbar.*(z(k+1)^3-z(k)^3);

                NT= NT + dT*Qbar*alphax.*(z(k+1) - z(k));
                NH= NH + dC*Qbar*betax.*(z(k+1) - z(k));
                MT= MT + (dT/2)*Qbar*alphax.*(z(k+1)^2-z(k)^2);
                MH= MH + (dC/2)*Qbar*betax.*(z(k+1)^2-z(k)^2);
            end
        end


        F=[A B;B D];

        T=[NT*1.0e3;MT*1.0e3];
        H=[NH*1.0e3;MH*1.0e3];

        epsilonF = inv(F) * N ;
        KF = [epsilonF(4:6,1)];
        epsilon0F = [epsilonF(1:3,1)];

        epsilonT= (inv(F) * T)*1.0e-3;
        KT = [epsilonT(4:6,1)];
        epsilon0T = [epsilonT(1:3,1)];


        epsilonH= (inv(F) * H)*1.0e-3;
        KH = [epsilonH(4:6,1)];
        epsilon0H = [epsilonH(1:3,1)];

        Sigma_F = zeros(3, n_lamina);
        Sigma_RT= zeros(3,n_lamina);
        Sigma_T = zeros(3, n_lamina);
        Sigma_RH= zeros(3,n_lamina);
        Sigma_H = zeros(3, n_lamina);


        if even==1
            for j=1:n_lamina
               
                     T=  [cosd(theta(j))^2, sind(theta(j))^2, 2*sind(theta(j))*cosd(theta(j));
                        sind(theta(j))^2, cosd(theta(j))^2, -2*sind(theta(j))*cosd(theta(j));
                        -sind(theta(j))*cosd(theta(j)), sind(theta(j))*cosd(theta(j)), cosd(theta(j))^2 - sind(theta(j))^2];

                    alpha_x= inv(T)*alpha1;
                    alphax=[alpha_x(1,1);alpha_x(2,1);2*alpha_x(3,1)];
                    beta_x= inv(T)*beta1;
                    betax=[beta_x(1,1);beta_x(2,1);2*beta_x(3,1)];

                    Qbar_11= Q11*(cosd(theta(j)))^4 + 2*(Q12+2*Q66)*(cosd(theta(j))^2) * (sind(theta(j)))^2 + Q22*(sind(theta(j)))^4;
                    Qbar_12= (Q11+Q22-4*Q66)*(cosd(theta(j)))^2 * (sind(theta(j)))^2 + Q12*((cosd(theta(j)))^4 + (sind(theta(j)))^4);
                    Qbar_16= (Q11-Q12-2*Q66)*cosd(theta(j))^3 * sind(theta(j)) - (Q22-Q12-2*Q66)*cosd(theta(j))*sind(theta(j))^3;
                    Qbar_21= (Q11+Q22-4*Q66)*cosd(theta(j))^2 * sind(theta(j))^2 + Q12*(cosd(theta(j))^4 + sind(theta(j))^4);
                    Qbar_22= Q11*sind(theta(j))^4 + 2*(Q12+2*Q66)*cosd(theta(j))^2 * sind(theta(j))^2 + Q22*cosd(theta(j))^4;
                    Qbar_26= (Q11-Q12-2*Q66)*sind(theta(j))^3 * cosd(theta(j))-(Q22-Q12-2*Q66)*cosd(theta(j))^3 * sind(theta(j));
                    Qbar_61= (Q11-Q12-2*Q66)*cosd(theta(j))^3 * sind(theta(j)) - (Q22-Q12-2*Q66)*cosd(theta(j))*sind(theta(j))^3;
                    Qbar_62= (Q11-Q12-2*Q66)*sind(theta(j))^3*cosd(theta(j))-(Q22-Q12-2*Q66)*cosd(theta(j))^3*sind(theta(j));
                    Qbar_66= (Q11+Q22-2*Q12-2*Q66)*cosd(theta(j))^2 * sind(theta(j))^2 + Q66*(cosd(theta(j))^4+sind(theta(j))^4);

                    Qbar =[Qbar_11 Qbar_12 Qbar_16; Qbar_21 Qbar_22 Qbar_26 ; Qbar_61 Qbar_62 Qbar_66];
 if j<=(n_lamina/2)
                    epsilon_xF = epsilon0F + z(j)*KF;
                    epsilon_1F= T*epsilon_xF ;
                    Sigma_k = Q * epsilon_1F;
                    % Accumulate the contribution
                    Sigma_F(:,j)=Sigma_F(:,j)+Sigma_k;

                    epsilon_x = epsilon0T + z(j)*KT;
                    epsilon_xT= dT*alphax;
                    epsilon_R = (epsilon_x) - (epsilon_xT);
                    Sigma_xR = Qbar * epsilon_R*1.0e3;
                    Sigma_1R =T*Sigma_xR;
                    % Accumulate the contribution
                    Sigma_RT(:,j) = Sigma_RT(:,j) + Sigma_1R;

                    epsilon_xh = epsilon0H + z(j)*KH;
                    epsilon_xH= dC*betax;
                    epsilon_RH = (epsilon_xh) - (epsilon_xH);
                    Sigma_xRH = Qbar * epsilon_RH*1.0e3;
                    Sigma_1RH =T*Sigma_xRH;
                    Sigma_RH(:,j) = Sigma_RH(:,j) + Sigma_1RH;
 elseif j>((n/2)+1)
     a=j+1;
             epsilon_xF = epsilon0F + z(a)*KF;
                    epsilon_1F= T*epsilon_xF ;
                    Sigma_k = Q * epsilon_1F;
                    % Accumulate the contribution
                    Sigma_F(:,j)=Sigma_F(:,j)+Sigma_k;

                    epsilon_x = epsilon0T + z(a)*KT;
                    epsilon_xT= dT*alphax;
                    epsilon_R = (epsilon_x) - (epsilon_xT);
                    Sigma_xR = Qbar * epsilon_R*1.0e3;
                    Sigma_1R =T*Sigma_xR;
                    % Accumulate the contribution
                    Sigma_RT(:,j) = Sigma_RT(:,j) + Sigma_1R;

                    epsilon_xh = epsilon0H + z(a)*KH;
                    epsilon_xH= dC*betax;
                    epsilon_RH = (epsilon_xh) - (epsilon_xH);
                    Sigma_xRH = Qbar * epsilon_RH*1.0e3;
                    Sigma_1RH =T*Sigma_xRH;
                    Sigma_RH(:,j) = Sigma_RH(:,j) + Sigma_1RH;
 end
            end
        else
            v=1;
            for j=1:n_lamina

                T=  [cosd(theta(j))^2, sind(theta(j))^2, 2*sind(theta(j))*cosd(theta(j));
                    sind(theta(j))^2, cosd(theta(j))^2, -2*sind(theta(j))*cosd(theta(j));
                    -sind(theta(j))*cosd(theta(j)), sind(theta(j))*cosd(theta(j)), cosd(theta(j))^2 - sind(theta(j))^2];

                alpha_x= inv(T)*alpha1;
                alphax=[alpha_x(1,1);alpha_x(2,1);2*alpha_x(3,1)];
                beta_x= inv(T)*beta1;
                betax=[beta_x(1,1);beta_x(2,1);2*beta_x(3,1)];

                Qbar_11= Q11*(cosd(theta(j)))^4 + 2*(Q12+2*Q66)*(cosd(theta(j))^2) * (sind(theta(j)))^2 + Q22*(sind(theta(j)))^4;
                Qbar_12= (Q11+Q22-4*Q66)*(cosd(theta(j)))^2 * (sind(theta(j)))^2 + Q12*((cosd(theta(j)))^4 + (sind(theta(j)))^4);
                Qbar_16= (Q11-Q12-2*Q66)*cosd(theta(j))^3 * sind(theta(j)) - (Q22-Q12-2*Q66)*cosd(theta(j))*sind(theta(j))^3;
                Qbar_21= (Q11+Q22-4*Q66)*cosd(theta(j))^2 * sind(theta(j))^2 + Q12*(cosd(theta(j))^4 + sind(theta(j))^4);
                Qbar_22= Q11*sind(theta(j))^4 + 2*(Q12+2*Q66)*cosd(theta(j))^2 * sind(theta(j))^2 + Q22*cosd(theta(j))^4;
                Qbar_26= (Q11-Q12-2*Q66)*sind(theta(j))^3 * cosd(theta(j))-(Q22-Q12-2*Q66)*cosd(theta(j))^3 * sind(theta(j));
                Qbar_61= (Q11-Q12-2*Q66)*cosd(theta(j))^3 * sind(theta(j)) - (Q22-Q12-2*Q66)*cosd(theta(j))*sind(theta(j))^3;
                Qbar_62= (Q11-Q12-2*Q66)*sind(theta(j))^3*cosd(theta(j))-(Q22-Q12-2*Q66)*cosd(theta(j))^3*sind(theta(j));
                Qbar_66= (Q11+Q22-2*Q12-2*Q66)*cosd(theta(j))^2 * sind(theta(j))^2 + Q66*(cosd(theta(j))^4+sind(theta(j))^4);

                Qbar =[Qbar_11 Qbar_12 Qbar_16; Qbar_21 Qbar_22 Qbar_26 ; Qbar_61 Qbar_62 Qbar_66];
                if j<(length(z)/2)
                    epsilon_xF = epsilon0F + z(v)*KF;
                    epsilon_1F= T*epsilon_xF ;
                    Sigma_k = Q * epsilon_1F;
                    % Accumulate the contribution
                    Sigma_F(:,j)=Sigma_F(:,j)+Sigma_k;

                    epsilon_x = epsilon0T + z(v)*KT;
                    epsilon_xT= dT*alphax;
                    epsilon_R = (epsilon_x) - (epsilon_xT);
                    Sigma_xR = Qbar * epsilon_R*1.0e3;
                    Sigma_1R =T*Sigma_xR;
                    % Accumulate the contribution
                    Sigma_RT(:,j) = Sigma_RT(:,j) + Sigma_1R;

                    epsilon_xh = epsilon0H + z(v)*KH;
                    epsilon_xH= dC*betax;
                    epsilon_RH = (epsilon_xh) - (epsilon_xH);
                    Sigma_xRH = Qbar * epsilon_RH*1.0e3;
                    Sigma_1RH =T*Sigma_xRH;
                    Sigma_RH(:,j) = Sigma_RH(:,j) + Sigma_1RH;
                    v=v+1;
                elseif j>=((length(z)/2)+1)
                     

                    epsilon_xF = epsilon0F + z(v)*KF;
                    epsilon_1F= T*epsilon_xF ;
                    Sigma_k = Q * epsilon_1F;
                    % Accumulate the contribution
                    Sigma_F(:,j)=Sigma_F(:,j)+Sigma_k;

                    epsilon_x = epsilon0T + z(v)*KT;
                    epsilon_xT= dT*alphax;
                    epsilon_R = (epsilon_x) - (epsilon_xT);
                    Sigma_xR = Qbar * epsilon_R*1.0e3;
                    Sigma_1R =T*Sigma_xR;
                    % Accumulate the contribution
                    Sigma_RT(:,j) = Sigma_RT(:,j) + Sigma_1R;

                    epsilon_xh = epsilon0H + z(v)*KH;
                    epsilon_xH= dC*betax;
                    epsilon_RH = (epsilon_xh) - (epsilon_xH);
                    Sigma_xRH = Qbar * epsilon_RH*1.0e3;
                    Sigma_1RH =T*Sigma_xRH;
                    Sigma_RH(:,j) = Sigma_RH(:,j) + Sigma_1RH;

                elseif  j==(length(z)/2)
                
                    for i=1:2
                        epsilon_xF = epsilon0F + z(v)*KF;
                        epsilon_1F= T*epsilon_xF ;
                        Sigma_k = Q * epsilon_1F;
                        % Accumulate the contribution
                        Sigma_F(:,j)=Sigma_F(:,j)+Sigma_k;

                        epsilon_x = epsilon0T + z(v)*KT;
                        epsilon_xT= dT*alphax;
                        epsilon_R = (epsilon_x) - (epsilon_xT);
                        Sigma_xR = Qbar * epsilon_R*1.0e3;
                        Sigma_1R =T*Sigma_xR;
                        % Accumulate the contribution
                        Sigma_RT(:,j) = Sigma_RT(:,j) + Sigma_1R;

                        epsilon_xh = epsilon0H + z(v)*KH;
                        epsilon_xH= dC*betax;
                        epsilon_RH = (epsilon_xh) - (epsilon_xH);
                        Sigma_xRH = Qbar * epsilon_RH*1.0e3;
                        Sigma_1RH =T*Sigma_xRH;
                        Sigma_RH(:,j) = Sigma_RH(:,j) + Sigma_1RH;
                        v=v+1;

                    end
                end

            end
        end

        for h=1:n_lamina
            if any(plyfail==theta(h))
                % continue
                sigma(:,h)=zeros(3,1);
            else

                if Sigma_F(1,h)>0
                    sigma(1,h)= Sigma_F(1,h)/sigma1_t;
                else
                    sigma(1,h)=Sigma_F(1,h)/sigma1_c;
                end
                if Sigma_F(2,h)>0
                    sigma(2,h)=Sigma_F(2,h)/sigma2_t;
                else
                    sigma(2,h)=Sigma_F(2,h)/sigma2_c;
                end
                sigma(3,h)= Sigma_F(3,h)/tau12;
            end
        end

        % disp(Sigma_F);
        % disp(Sigma_R);
        disp(sigma);
        % Find the maximum value and its column index from the entire sigma matrix
        [max_ratio, max_column_index] = max(abs(sigma(:)));

        % Convert linear index to row and column indices
        [row_index, column_index] = ind2sub(size(sigma), max_column_index);

        fprintf('Lamina Failed: %d \n',column_index);
        fprintf(' %d degree ply fails in %s  since max_ratio is for (%d,%d)\n',theta(column_index),char(failure_criteria(row_index)),row_index,column_index);

        % Display the maximum value and its column index
        fprintf('Maximum ratio: %f\n', max_ratio);
        % fprintf('max stress ratio comes out for %d degree ply \n', theta(column_index));

        for l=1:n_lamina
            if l==row_index
                Total_stress = Sigma_F + Sigma_RT + Sigma_RH;
                % Initialize strength
                strength=0;
                % Check each component of Total_stress
                for l=1:n_lamina
                    if sigma(row_index,column_index) > 0 && row_index==1
                        strength = sigma1_t;
                    elseif sigma(row_index,column_index) < 0 && row_index==1
                        strength =  sigma1_c;
                    end

                    if sigma(row_index,column_index) > 0 && row_index==2
                        strength =  sigma2_t;
                    elseif sigma(row_index,column_index) < 0 && row_index==2
                        strength =  sigma2_c;
                    end

                    if  sigma(row_index,column_index) > 0 && sigma(row_index,column_index) < 0 && row_index==3
                        strength= tau12;
                    end

                    stress_dueto_Nx= strength - Sigma_RT(row_index,column_index)-Sigma_RH(row_index,column_index);
                    FPF_LOAD = (N/(Sigma_F(row_index,column_index))*stress_dueto_Nx);
                end
            end
        end

        fprintf('failure load for %d degree ply is: %0.8f\n',theta(column_index),FPF_LOAD(1,1));


        plyfail = [plyfail, theta(column_index)];
        % plyfail(:,k_outer)= plyfail(:,k_outer)* theta(column_index);

        failure_LOAD(:,k_outer) = failure_LOAD(:,k_outer) + FPF_LOAD;

        % end
        N=FPF_LOAD;
    end

    disp('=============================================================================================')
    fprintf('FPF load is: %0.8f\n',failure_LOAD(1,1));
    fprintf('LPF load is: %0.8f\n',failure_LOAD(1,n));
    disp('======================================================================================')

