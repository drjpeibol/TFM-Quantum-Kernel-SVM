%------------------------------------------------------------------------
%  Seven_qbit_sample_generation_v0 : program to check the data generation 
%                                    explaind in  the paper J.R. Click 
%                                    "Covariant quantum kernels 
%                                     for data with group structure" 2022
%---------------------------------------------------------------------------

function Seven_qbit_sample_generation_v0
%Symbolic calculations
symbolic = 0;  %Use 1 to get the body of the U function
test_number = 1; % Loops to test functions
samples = 20; % Number of samples to generate
eps = 0.01; %Variance of the random error applyed
seed = 1111;

rng(seed);

% Identity 
I = [ 1 0 
      0 1 ];
% Pauli Matrices
X = [ 0 1
      1 0 ];
Y = [ 0  -1i
      1i   0];
Z = [  1  0 
       0 -1 ];


%-------------------------------------------------------------------------
% Symbolic generation of
% kron( D(theta_1, theta_2,0 ), D(theta_3, theta_4,0))

    if symbolic == 1
    syms theta_1
    Rot1 = [ cos(theta_1/2)    -1i*sin(theta_1/2);
            -1i*sin(theta_1/2)    cos(theta_1/2)  ];

    syms theta_2
    Rot2 = [ exp(-1i*theta_2/2)    0
            0              exp(1i*theta_2/2) ]; 
    syms theta_3
    Rot3 = [ cos(theta_3/2)    -1i*sin(theta_3/2);
            -1i*sin(theta_3/2)    cos(theta_3/2)  ];

    syms theta_4
    Rot4 = [ exp(-1i*theta_4/2)    0
            0              exp(1i*theta_4/2) ];

    Result = kron(Rot1*Rot2, Rot3*Rot4);

    fprintf ('Tensor product decompositon:\n');

    fprintf("M(1,1) = %s; \nM(1,2)= %s;\nM(1,3) = %s;\nM(1,4) = %s;\n",....
             Result(1,1), Result(1,2), Result(1,3), Result(1,4)); 
    fprintf("M(2,1) = %s; \nM(2,2)= %s;\nM(2,3) = %s;\nM(2,4) = %s;\n",....
             Result(2,1), Result(2,2), Result(2,3), Result(2,4));
    fprintf("M(3,1) = %s; \nM(3,2)= %s;\nM(3,3) = %s;\nM(3,4) = %s;\n",....
             Result(3,1), Result(3,2), Result(3,3), Result(3,4)); 
    fprintf("M(4,1) = %s; \nM(4,2)= %s;\nM(4,3) = %s;\nM(4,4) = %s;\n",....
             Result(4,1), Result(4,2), Result(4,3), Result(4,4));    
                 
    end
%-------------------------------------------------------------------------

%-------------------------------------------------------------------- 
%
%  Use of Five qbits for the architecture:
%          0 --- 1 --- 2 
%                |
%                3
%                |
%          4 --- 5 --- 6
%
%  Therefore, we have six cases:
%  X0Z1
%  X1Z2
%  X1Z3
%  X3Z5
%  X4Z5
%  X5Z6
%--------------------------------------------------------------------

coset = [ 'C-' 
          'C+' ];

for k = 1:test_number
 for j = 1:2
    
    % Indicate the coset
    fprintf ('Coset %s\n', coset(j,:));
    
    % Generate random angles  
    theta_1 = random('uniform', -pi/2, pi/2);
    theta_2 = random('uniform', -pi/2, pi/2);
    theta_3 = random('uniform', -pi/2, pi/2);
    theta_4 = random('uniform', -pi/2, pi/2);
    theta_5 = random('uniform', -pi/2, pi/2);
    theta_6 = random('uniform', -pi/2, pi/2);
    theta_7 = random('uniform', -pi/2, pi/2);
    theta_8 = random('uniform', -pi/2, pi/2);
    theta_9 = random('uniform', -pi/2, pi/2);
    theta_10 = random('uniform', -pi/2, pi/2);
    theta_11 = random('uniform', -pi/2, pi/2);
    theta_12 = random('uniform', -pi/2, pi/2);
    theta_13 = random('uniform', -pi/2, pi/2);
    theta_14 = random('uniform', -pi/2, pi/2);
    
    % Identity 
    Sample = [ theta_1, theta_2, ...
               theta_3, theta_4, ...
               theta_5, theta_6, ...
               theta_7, theta_8, ...
               theta_9, theta_10,...
               theta_11,theta_12,...
               theta_13,theta_14];
    fprintf ('Sample R*I    =  %.2f, %.2f, %.2f, %.2f, %.2f,%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n',...
              Sample); 
    %----------------------------------------------------------------
    % First case: X0Z1
    % ---  X --- D(0, theta_1, theta_2, 0) --- 
    % ---  Z --- D(0, theta_3, theta_4, 0) --- 
    % ---- I ----D(0, theta_5, theta_6, 0) ---
    % ---- I ----D(0, theta_7, theta_8, 0) ---
    % ---- I ----D(0, theta_9, theta_10, 0) ---
    % ---- I ----D(0, theta_11, theta_12, 0) ---
    % ---- I ----D(0, theta_13, theta_14, 0) ---

    
    % Rotations  
    R = kron(D([0, theta_1, theta_2, 0], I, X, Z, X), ...
         kron(D([0, theta_3, theta_4, 0], I, X, Z, X), ...
           kron( D([0, theta_5, theta_6, 0], I, X, Z, X), ...
             kron( D([0, theta_7, theta_8, 0], I, X, Z, X), ...
                kron( D([0, theta_9, theta_10, 0], I, X, Z, X), ... 
                   kron( D([0, theta_11, theta_12, 0], I, X, Z, X),   ...
                         D([0, theta_13, theta_14, 0], I, X, Z, X) ...
                                                             ))))));
    % Stabilizer 
    X0Z1 = kron(X, ...
            kron(Z, ...
                kron(I, ...
                    kron(I, ...
                        kron(I, ...
                            kron (I, ...
                                    I...
                    ))))));
    
    A = R*X0Z1;
    
    % We parametrise the matrix by using first two qbits 
    M0 = kron(D([0, theta_1, theta_2, 0], I, X, Z, X), ...
           D([0, theta_3, theta_4, 0], I, X, Z, X)     );
    
    S0= kron(X,Z);
    B0 = M0*S0;
    
    % Recalculate angles 
    estimated_angles = DxD_Euler_angles (B0);
    

    % Equivalent Matrix 
    % --- D(0, estimated_angles(1),estimated_angles(2), 0) --- 
    % --- D(0, estimated_angles(3),estimated_angles(4), 0) --- 
    % ----D(0, theta_5, theta_6, 0) --------------------------
    % ----D(0, theta_7, theta_8, 0) --------------------------
    % ----D(0, theta_9, theta_10, 0) -------------------------
    % ----D(0, theta_11, theta_12, 0) ------------------------
    % ----D(0, theta_13, theta_14, 0) ------------------------
    
                                                   
    A0 = kron(D([0,estimated_angles(1),estimated_angles(2), 0], I,X,Z,X), ...
          kron( D([0,estimated_angles(3),estimated_angles(4), 0], I,X,Z,X), ...
           kron( D([0, theta_5, theta_6, 0], I, X, Z, X), ...
             kron( D([0, theta_7, theta_8, 0], I, X, Z, X), ...
                kron( D([0, theta_9, theta_10, 0], I, X, Z, X), ... 
                   kron( D([0, theta_11, theta_12, 0], I, X, Z, X),   ...
                         D([0, theta_13, theta_14, 0], I, X, Z, X) ...
                                                             ))))));                                            
    
    %test 
    if norm (A-A0)  >= 0.001
         fprintf ('Matrix angles estgimation warning\n');
         return
    else

    end
    
    % Therefore, sample is 
    Sample0 = [ estimated_angles(1),estimated_angles(2), ...
               estimated_angles(3),estimated_angles(4), ...
               theta_5, theta_6, ...
               theta_7, theta_8, ...
               theta_9, theta_10,...
               theta_11,theta_12,...
               theta_13,theta_14];
    fprintf ('Sample R*X0Z1 =  %.2f, %.2f, %.2f, %.2f, %.2f,%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n',...
              Sample0);
    if j == 1
        writematrix ([Sample0,j],'Data_7qubit.csv','WriteMode','overwrite');
    else 
        writematrix ([Sample0,j],'Data_7qubit.csv','WriteMode','append');
    end 
    %----
    %  End of first case 
    %
    %------------------------------------------------------------------
    % Second case: X1Z2 
    % ---  I --- D(0, theta_1, theta_2, 0) --- 
    % ---  X --- D(0, theta_3, theta_4, 0) --- 
    % ---- Z ----D(0, theta_5, theta_6, 0) ---
    % ---- I ----D(0, theta_7, theta_8, 0) ---
    % ---- I ----D(0, theta_9, theta_10, 0) ---
    % ---- I ----D(0, theta_11, theta_12, 0) ---
    % ---- I ----D(0, theta_13, theta_14, 0) ---
    
    % Stabilizer 
    X1Z2 = kron(I, ...
            kron(X, ...
                kron(Z, ...
                    kron(I, ...
                        kron(I, ...
                            kron (I, ...
                                    I...
                    ))))));
                  
    A = R*X1Z2;
   
    % Parametrise q1 and q2
    M1 = kron(D([0, theta_3, theta_4, 0], I, X, Z, X), ...
           D([0, theta_5, theta_6, 0], I, X, Z, X)     );
    S1 = kron(X,Z);
    B1 = M1*S1;
    % Recalculate angles 
    estimated_angles = DxD_Euler_angles (B1);
   
    % Equivalent Matrix 
    % ----D(0, theta_1, theta_2, 0) --------------------------
    % --- D(0, estimated_angles(1),estimated_angles(2), 0) --- 
    % --- D(0, estimated_angles(3),estimated_angles(4), 0) --- 
    % ----D(0, theta_7, theta_8, 0) --------------------------
    % ----D(0, theta_9, theta_10, 0) -------------------------
    % ---- I ----D(0, theta_11, theta_12, 0) -----------------
    % ---- I ----D(0, theta_13, theta_14, 0) -----------------

    A1 = kron( D([0, theta_1, theta_2, 0], I, X, Z, X), ...  
          kron(D([0,estimated_angles(1),estimated_angles(2), 0], I,X,Z,X), ...
           kron( D([0,estimated_angles(3),estimated_angles(4), 0], I,X,Z,X), ...
             kron( D([0, theta_7, theta_8, 0], I, X, Z, X), ...
                kron( D([0, theta_9, theta_10, 0], I, X, Z, X), ... 
                   kron( D([0, theta_11, theta_12, 0], I, X, Z, X),   ...
                         D([0, theta_13, theta_14, 0], I, X, Z, X) ...
                                                             ))))));                                              
    
    %test 
    if norm (A-A1)  >= 0.001
         fprintf ('Matrix angles estgimation warning\n');
         return
    end
        
    % Therefore, sample is 
    Sample_1 = [ theta_1, theta_2,...
               estimated_angles(1),estimated_angles(2), ...
               estimated_angles(3),estimated_angles(4), ...
               theta_7, theta_8, ...
               theta_9, theta_10,...
               theta_11,theta_12,...
               theta_13,theta_14];
    fprintf ('Sample R*X1Z2 =  %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n',...
              Sample_1);

    writematrix ([Sample_1,j],'Data_7qubit.csv','WriteMode','append');

    %
    % End of second case
    %-------------------------------------------------------------------
    % Third case: X1Z3
    % ---  I --- D(0, theta_1, theta_2, 0) --- 
    % ---  X --- D(0, theta_3, theta_4, 0) --- 
    % ---- I ----D(0, theta_5, theta_6, 0) ---
    % ---- Z ----D(0, theta_7, theta_8, 0) ---
    % ---- I ----D(0, theta_9, theta_10, 0) ---
    % ---- I ----D(0, theta_11, theta_12, 0) ---
    % ---- I ----D(0, theta_13, theta_14, 0) ---
    
    % Stabilizer 
    X1Z3 = kron(I, ...
            kron(X, ...
                kron(I, ...
                    kron(Z, ...
                        kron(I, ...
                            kron (I, ...
                                    I...
                    ))))));
                  
    A = R*X1Z3;
    
    % Interchange q2 and q3 
    % ---  I --- D(0, theta_1, theta_2, 0) --- 
    % ---  X --- D(0, theta_3, theta_4, 0) --- 
    % ---- Z ----D(0, theta_7, theta_8, 0) ---
    % ---- I ----D(0, theta_5, theta_6, 0) ---
    % ---- I ----D(0, theta_9, theta_10, 0) ---
    % ---- I ----D(0, theta_11, theta_12, 0) ---
    % ---- I ----D(0, theta_13, theta_14, 0) ---

    P = [ 1 0 0 0 
         0 0 1 0  
         0 1 0 0  
         0 0 0 1   ]; % Interchange two qbits

     P = kron(I, ...
           kron(I, ...
             kron(P, ...
                kron(I, ...
                    kron(I, ...
                            I...
                         )))));
    
        
    
    % Parametrise q1 and q3
    M1 = kron(D([0, theta_3, theta_4, 0], I, X, Z, X), ...
           D([0, theta_7, theta_8, 0], I, X, Z, X)     );

    S1 = kron(X,Z);
    B1 = M1*S1;
    % Recalculate angles 
    estimated_angles = DxD_Euler_angles (B1);
    
    % Equivalent Matrix 
    % ----D(0, theta_1, theta_2, 0) --------------------------
    % --- D(0, estimated_angles(1),estimated_angles(2), 0) --- 
    % --- D(0, estimated_angles(3),estimated_angles(4), 0) --- 
    % ----D(0, theta_5, theta_6, 0) --------------------------
    % ----D(0, theta_9, theta_10, 0) -------------------------
    % ----D(0, theta_10, theta_11, 0) --------------------------
    % ----D(0, theta_12, theta_13, 0) -------------------------

    A1 = kron( D([0, theta_1, theta_2, 0], I, X, Z, X), ...  
          kron(D([0,estimated_angles(1),estimated_angles(2), 0], I,X,Z,X), ...
           kron( D([0,estimated_angles(3),estimated_angles(4), 0], I,X,Z,X), ...
             kron( D([0, theta_5, theta_6, 0], I, X, Z, X), ... 
                kron( D([0, theta_9, theta_10, 0], I, X, Z, X), ... 
                   kron( D([0, theta_11, theta_12, 0], I, X, Z, X),   ...
                         D([0, theta_13, theta_14, 0], I, X, Z, X) ...
                                                             ))))));  
     
    % undo swap 
    A1=P*A1*P;

    
    %test 
    if norm (A-A1)  >= 0.001
         fprintf ('Matrix angles estimation warning\n');
         return;
    end
        
    % Therefore, sample is 
    Sample_2 = [ theta_1, theta_2,...
               estimated_angles(1),estimated_angles(2), ...
               theta_5, theta_6, ...
               estimated_angles(3),estimated_angles(4), ...
               theta_9, theta_10,...
               theta_11,theta_12,...
               theta_13,theta_14];
    fprintf ('Sample R*X1Z3 =  %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n',...
              Sample_2);

    writematrix ([Sample_2,j],'Data_7qubit.csv','WriteMode','append');

    % End of third case 
    %---------------------------------------------------------------------
    % Fourth case: X3Z5 
    % ---  I --- D(0, theta_1, theta_2, 0) --- 
    % ---  I --- D(0, theta_3, theta_4, 0) --- 
    % ---- I ----D(0, theta_5, theta_6, 0) ---
    % ---- X ----D(0, theta_7, theta_8, 0) ---
    % ---- I ----D(0, theta_9, theta_10, 0) ---
    % ---- Z ----D(0, theta_11, theta_12, 0) ---
    % ---- I ----D(0, theta_13, theta_14, 0) ---
    
 
    % Stabilizer 
    X3Z5 = kron(I, ...
            kron(I, ...
                kron(I, ...
                    kron(X, ...
                        kron(I, ...
                            kron(Z, ...
                                    I...
                    ))))));
                  
    A = R*X3Z5;
    
    % Interchange q4 and q5 
    % ---  I --- D(0, theta_1, theta_2, 0) --- 
    % ---  I --- D(0, theta_3, theta_4, 0) --- 
    % ---- I ----D(0, theta_5, theta_6, 0) ---
    % ---- X ----D(0, theta_7, theta_8, 0) ---
    % ---- Z ----D(0, theta_11, theta_12, 0) ---
    % ---- I ----D(0, theta_9, theta_10, 0) ---
    % ---- I ----D(0, theta_13, theta_14, 0) ---

    
    P = [ 1 0 0 0 
         0 0 1 0  
         0 1 0 0  
         0 0 0 1   ]; % Interchange two qbits

     P = kron(I, ...
           kron(I, ...
             kron(I, ...
                kron(I, ...
                    kron(P, ...
                            I...
                         )))));
    
        
    
    % Parametrise q4 and q5
    M1 = kron(D([0, theta_7, theta_8, 0], I, X, Z, X), ...
           D([0, theta_11, theta_12, 0], I, X, Z, X)     );

    S1 = kron(X,Z);
    B1 = M1*S1;
    % Recalculate angles 
    estimated_angles = DxD_Euler_angles (B1);
    

    % Equivalent Matrix 
    % ----D(0, theta_1, theta_2, 0) --------------------------
    % ----D(0, theta_3, theta_4, 0) --------------------------
    % ----D(0, theta_5, theta_6, 0) -------------------------
    % --- D(0, estimated_angles(1),estimated_angles(2), 0) --- 
    % --- D(0, estimated_angles(3),estimated_angles(4), 0) --- 
    % ----D(0, theta_9, theta_10, 0) --------------------------
    % ----D(0, theta_12, theta_13, 0) -------------------------

    A1 = kron( D([0, theta_1, theta_2, 0], I, X, Z, X), ...  
          kron(D([0,theta_3,theta_4, 0], I, X, Z, X), ...
           kron( D([0,theta_5,theta_6, 0], I,X,Z,X), ...
             kron( D([0, estimated_angles(1), estimated_angles(2), 0], I, X, Z, X), ... 
                kron( D([0, estimated_angles(3), estimated_angles(4), 0], I, X, Z, X), ... 
                   kron( D([0, theta_9, theta_10, 0], I, X, Z, X), ...
                         D([0, theta_13, theta_14, 0], I, X, Z, X)...
                                                             ))))));  
     
    % undo swap 
    A1=P*A1*P;
    
    %test 
    if norm (A-A1)  >= 0.001
         fprintf ('Matrix X3Z5 angles estgimation warning\n');
         return;
    end
        
    % Therefore, sample is 
    Sample_3 = [ theta_1, theta_2,...
               theta_3,theta_4, ...
               theta_5, theta_6, ...
               estimated_angles(1),estimated_angles(2), ...
               theta_9, theta_10,...
               estimated_angles(3),estimated_angles(4),...
               theta_13,theta_14];
    fprintf ('Sample R*X3Z5 =  %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n',...
              Sample_3);

    writematrix ([Sample_3,j],'Data_7qubit.csv','WriteMode','append');


    % End of fourth case 
    %---------------------------------------------------------------------
    % Fifth case: X4Z5 
    % ---  I --- D(0, theta_1, theta_2, 0) --- 
    % ---  I --- D(0, theta_3, theta_4, 0) --- 
    % ---- I ----D(0, theta_5, theta_6, 0) ---
    % ---- I ----D(0, theta_7, theta_8, 0) ---
    % ---- X ----D(0, theta_9, theta_10, 0) ---
    % ---- Z ----D(0, theta_11, theta_12, 0) ---
    % ---- I ----D(0, theta_13, theta_14, 0) ---
    
    % Stabilizer 
    X4Z5 = kron(I, ...
            kron(I, ...
                kron(I, ...
                    kron(I, ...
                        kron(X, ...
                            kron (Z, ...
                                    I...
                    ))))));
                  
    A = R*X4Z5;
   
    % Parametrise q1 and q2
    M1 = kron(D([0, theta_9, theta_10, 0], I, X, Z, X), ...
           D([0, theta_11, theta_12, 0], I, X, Z, X)     );
    S1 = kron(X,Z);
    B1 = M1*S1;
    % Recalculate angles 
    estimated_angles = DxD_Euler_angles (B1);
   
    % Equivalent Matrix 
    % ----D(0, theta_1, theta_2, 0) --------------------------
    % ----D(0, theta_3, theta_4, 0) --------------------------
    % ----D(0, theta_6, theta_5, 0) --------------------------
    % ----D(0, theta_7, theta_8, 0) --------------------------
    % --- D(0, estimated_angles(1),estimated_angles(2), 0) --- 
    % --- D(0, estimated_angles(3),estimated_angles(4), 0) --- 
    % ---- I ----D(0, theta_13, theta_14, 0) -----------------

    A1 = kron( D([0, theta_1, theta_2, 0], I, X, Z, X), ...  
          kron(D([0,theta_3,theta_4, 0], I,X,Z,X), ...
           kron( D([0,theta_5,theta_6, 0], I,X,Z,X), ...
             kron( D([0, theta_7, theta_8, 0], I, X, Z, X), ...
                kron( D([0, estimated_angles(1), estimated_angles(2), 0], I, X, Z, X), ... 
                   kron( D([0, estimated_angles(3), estimated_angles(4), 0], I, X, Z, X),   ...
                         D([0, theta_13, theta_14, 0], I, X, Z, X) ...
                                                             ))))));                                              
    
    %test 
    if norm (A-A1)  >= 0.001
         fprintf ('Matrix angles estgimation warning\n');
         return
    end
        
    % Therefore, sample is 
    Sample_4 = [ theta_1, theta_2,...
               theta_3,theta_4, ...
               theta_5,theta_6, ...
               theta_7, theta_8, ...
               estimated_angles(1),estimated_angles(2), ...
               estimated_angles(3),estimated_angles(4), ...
               theta_13,theta_14];
    fprintf ('Sample R*X4Z5 =  %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n',...
              Sample_4);

    writematrix ([Sample_4,j],'Data_7qubit.csv','WriteMode','append');

    % End of fifth case 
    %---------------------------------------------------------------------
    % Sixth case: X5Z6 
    % ---  I --- D(0, theta_1, theta_2, 0) --- 
    % ---  I --- D(0, theta_3, theta_4, 0) --- 
    % ---- I ----D(0, theta_5, theta_6, 0) ---
    % ---- I ----D(0, theta_7, theta_8, 0) ---
    % ---- I ----D(0, theta_9, theta_10, 0) ---
    % ---- X ----D(0, theta_11, theta_12, 0) ---
    % ---- Z ----D(0, theta_13, theta_14, 0) ---
    
    % Stabilizer 
    X5Z6 = kron(I, ...
            kron(I, ...
                kron(I, ...
                    kron(I, ...
                        kron(I, ...
                            kron (X, ...
                                    Z...
                    ))))));
                  
    A = R*X5Z6;
   
    % Parametrise q5 and q6
    M1 = kron(D([0, theta_11, theta_12, 0], I, X, Z, X), ...
           D([0, theta_13, theta_14, 0], I, X, Z, X)     );
    S1 = kron(X,Z);
    B1 = M1*S1;
    % Recalculate angles 
    estimated_angles = DxD_Euler_angles (B1);
   
    % Equivalent Matrix 
    % ----D(0, theta_1, theta_2, 0) --------------------------
    % ----D(0, theta_3, theta_4, 0) --------------------------
    % ----D(0, theta_6, theta_5, 0) --------------------------
    % ----D(0, theta_7, theta_8, 0) --------------------------
    % --- D(0, estimated_angles(1),estimated_angles(2), 0) --- 
    % --- D(0, estimated_angles(3),estimated_angles(4), 0) --- 
    % ---- I ----D(0, theta_13, theta_14, 0) -----------------

    A1 = kron( D([0, theta_1, theta_2, 0], I, X, Z, X), ...  
          kron(D([0,theta_3,theta_4, 0], I,X,Z,X), ...
           kron( D([0,theta_5,theta_6, 0], I,X,Z,X), ...
             kron( D([0, theta_7, theta_8, 0], I, X, Z, X), ...
                kron( D([0, theta_9, theta_10, 0], I, X, Z, X), ... 
                   kron( D([0, estimated_angles(1), estimated_angles(2), 0], I, X, Z, X),   ...
                         D([0, estimated_angles(3), estimated_angles(4), 0], I, X, Z, X) ...
                                                             ))))));                                              
    
    %test 
    if norm (A-A1)  >= 0.001
         fprintf ('Matrix angles estgimation warning\n');
         return
    end
        
    % Therefore, sample is 
    Sample_5 = [ theta_1, theta_2,...
               theta_3,theta_4, ...
               theta_5,theta_6, ...
               theta_7, theta_8, ...
               theta_9, theta_10, ...
               estimated_angles(1),estimated_angles(2), ...
               estimated_angles(3),estimated_angles(4)];
    fprintf ('Sample R*X5Z6 =  %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n',...
              Sample_5);

    writematrix ([Sample_5,j],'Data_7qubit.csv','WriteMode','append');
%---------------------------------------------------------------------
% Generate Dataset
    for r = 1:samples

%c_minus*s_1
        writematrix ([[Sample0(1),Sample0(2),Sample0(3),Sample0(4),Sample0(5),Sample0(6),Sample0(7),Sample0(8),Sample0(9),Sample0(10),Sample0(11),Sample0(12),Sample0(13),Sample0(14)] + ...
            [Sample0(1)*random("Normal",0,sqrt(eps)),Sample0(2)*random("Normal",0,sqrt(eps)),Sample0(3)*random("Normal",0,sqrt(eps)), ...
            Sample0(4)*random("Normal",0,sqrt(eps)),Sample0(5)*random("Normal",0,sqrt(eps)),Sample0(6)*random("Normal",0,sqrt(eps)), ...
            Sample0(7)*random("Normal",0,sqrt(eps)),Sample0(8)*random("Normal",0,sqrt(eps)),Sample0(9)*random("Normal",0,sqrt(eps)),Sample0(10)*random("Normal",0,sqrt(eps)), ...
            Sample0(11)*random("Normal",0,sqrt(eps)),Sample0(12)*random("Normal",0,sqrt(eps)),Sample0(13)*random("Normal",0,sqrt(eps)),Sample0(14)*random("Normal",0,sqrt(eps))] ...
            ,j],'Data_7qubit.csv','WriteMode','append');

        writematrix ([[Sample_1(1),Sample_1(2),Sample_1(3),Sample_1(4),Sample_1(5),Sample_1(6),Sample_1(7),Sample_1(8),Sample_1(9),Sample_1(10),Sample_1(11),Sample_1(12),Sample_1(13),Sample_1(14)] + ...
            [Sample_1(1)*random("Normal",0,sqrt(eps)),Sample_1(2)*random("Normal",0,sqrt(eps)),Sample_1(3)*random("Normal",0,sqrt(eps)), ...
            Sample_1(4)*random("Normal",0,sqrt(eps)),Sample_1(5)*random("Normal",0,sqrt(eps)),Sample_1(6)*random("Normal",0,sqrt(eps)), ...
            Sample_1(7)*random("Normal",0,sqrt(eps)),Sample_1(8)*random("Normal",0,sqrt(eps)),Sample_1(9)*random("Normal",0,sqrt(eps)),Sample_1(10)*random("Normal",0,sqrt(eps)), ...
            Sample_1(11)*random("Normal",0,sqrt(eps)),Sample_1(12)*random("Normal",0,sqrt(eps)),Sample_1(13)*random("Normal",0,sqrt(eps)),Sample_1(14)*random("Normal",0,sqrt(eps))] ...
            ,j],'Data_7qubit.csv','WriteMode','append');

        writematrix ([[Sample_2(1),Sample_2(2),Sample_2(3),Sample_2(4),Sample_2(5),Sample_2(6),Sample_2(7),Sample_2(8),Sample_2(9),Sample_2(10),Sample_2(11),Sample_2(12),Sample_2(13),Sample_2(14)] + ...
            [Sample_2(1)*random("Normal",0,sqrt(eps)),Sample_2(2)*random("Normal",0,sqrt(eps)),Sample_2(3)*random("Normal",0,sqrt(eps)), ...
            Sample_2(4)*random("Normal",0,sqrt(eps)),Sample_2(5)*random("Normal",0,sqrt(eps)),Sample_2(6)*random("Normal",0,sqrt(eps)), ...
            Sample_2(7)*random("Normal",0,sqrt(eps)),Sample_2(8)*random("Normal",0,sqrt(eps)),Sample_2(9)*random("Normal",0,sqrt(eps)),Sample_2(10)*random("Normal",0,sqrt(eps)), ...
            Sample_2(11)*random("Normal",0,sqrt(eps)),Sample_2(12)*random("Normal",0,sqrt(eps)),Sample_2(13)*random("Normal",0,sqrt(eps)),Sample_2(14)*random("Normal",0,sqrt(eps))] ...
            ,j],'Data_7qubit.csv','WriteMode','append');

        writematrix ([[Sample_3(1),Sample_3(2),Sample_3(3),Sample_3(4),Sample_3(5),Sample_3(6),Sample_3(7),Sample_3(8),Sample_3(9),Sample_3(10),Sample_3(11),Sample_3(12),Sample_3(13),Sample_3(14)] + ...
            [Sample_3(1)*random("Normal",0,sqrt(eps)),Sample_3(2)*random("Normal",0,sqrt(eps)),Sample_3(3)*random("Normal",0,sqrt(eps)), ...
            Sample_3(4)*random("Normal",0,sqrt(eps)),Sample_3(5)*random("Normal",0,sqrt(eps)),Sample_3(6)*random("Normal",0,sqrt(eps)), ...
            Sample_3(7)*random("Normal",0,sqrt(eps)),Sample_3(8)*random("Normal",0,sqrt(eps)),Sample_3(9)*random("Normal",0,sqrt(eps)),Sample_3(10)*random("Normal",0,sqrt(eps)), ...
            Sample_3(11)*random("Normal",0,sqrt(eps)),Sample_3(12)*random("Normal",0,sqrt(eps)),Sample_3(13)*random("Normal",0,sqrt(eps)),Sample_3(14)*random("Normal",0,sqrt(eps))] ...
            ,j],'Data_7qubit.csv','WriteMode','append');

        writematrix ([[Sample_4(1),Sample_4(2),Sample_4(3),Sample_4(4),Sample_4(5),Sample_4(6),Sample_4(7),Sample_4(8),Sample_4(9),Sample_4(10),Sample_4(11),Sample_4(12),Sample_4(13),Sample_4(14)] + ...
            [Sample_4(1)*random("Normal",0,sqrt(eps)),Sample_4(2)*random("Normal",0,sqrt(eps)),Sample_4(3)*random("Normal",0,sqrt(eps)), ...
            Sample_4(4)*random("Normal",0,sqrt(eps)),Sample_4(5)*random("Normal",0,sqrt(eps)),Sample_4(6)*random("Normal",0,sqrt(eps)), ...
            Sample_4(7)*random("Normal",0,sqrt(eps)),Sample_4(8)*random("Normal",0,sqrt(eps)),Sample_4(9)*random("Normal",0,sqrt(eps)),Sample_4(10)*random("Normal",0,sqrt(eps)), ...
            Sample_4(11)*random("Normal",0,sqrt(eps)),Sample_4(12)*random("Normal",0,sqrt(eps)),Sample_4(13)*random("Normal",0,sqrt(eps)),Sample_4(14)*random("Normal",0,sqrt(eps))] ...
            ,j],'Data_7qubit.csv','WriteMode','append');

        writematrix ([[Sample_5(1),Sample_5(2),Sample_5(3),Sample_5(4),Sample_5(5),Sample_5(6),Sample_5(7),Sample_5(8),Sample_5(9),Sample_5(10),Sample_5(11),Sample_5(12),Sample_5(13),Sample_5(14)] + ...
            [Sample_5(1)*random("Normal",0,sqrt(eps)),Sample_5(2)*random("Normal",0,sqrt(eps)),Sample_5(3)*random("Normal",0,sqrt(eps)), ...
            Sample_5(4)*random("Normal",0,sqrt(eps)),Sample_5(5)*random("Normal",0,sqrt(eps)),Sample_5(6)*random("Normal",0,sqrt(eps)), ...
            Sample_5(7)*random("Normal",0,sqrt(eps)),Sample_5(8)*random("Normal",0,sqrt(eps)),Sample_5(9)*random("Normal",0,sqrt(eps)),Sample_5(10)*random("Normal",0,sqrt(eps)), ...
            Sample_5(11)*random("Normal",0,sqrt(eps)),Sample_5(12)*random("Normal",0,sqrt(eps)),Sample_5(13)*random("Normal",0,sqrt(eps)),Sample_5(14)*random("Normal",0,sqrt(eps))] ...
            ,j],'Data_7qubit.csv','WriteMode','append');

    end
 end  
end

end



%------------------------------------------------------------------------
% DxD_Euler_angles : Caluclates the Euler angles for a SU(4) = DxD matitx
%------------------------------------------------------------------------
function angles = DxD_Euler_angles (M)

%M is SU(4)
I= [ 1 0 0 0 
     0 1 0 0 
     0 0 1 0 
     0 0 0 1 ];
 if norm(M*M'-I) >= 0.001    % not U4
     angles=[0 0 0 0];
     return
 elseif det(M)-1 >= 0.001  % not SU(4)
     angles = [0 0 0 0];
     return
 end
         

% Using M explicit expressions 
% M(1,2)= -cos(theta_1/2)*sin(theta_3/2)*exp(-(theta_2*1i)/2)*exp((theta_4*1i)/2)*1i;
% M(1,1) = cos(theta_1/2)*cos(theta_3/2)*exp(-(theta_2*1i)/2)*exp(-(theta_4*1i)/2); 
theta_3 =   2*atan( abs(M(1,2)/M(1,1)) );

% M(1,4) = -sin(theta_1/2)*sin(theta_3/2)*exp((theta_2*1i)/2)*exp((theta_4*1i)/2);
% M(2,1) = -cos(theta_1/2)*sin(theta_3/2)*exp(-(theta_2*1i)/2)*exp(-(theta_4*1i)/2)*1i; 
theta_1 =   2*atan( abs(M(1,4)/M(2,1)) );

% M(2,2)= cos(theta_1/2)*cos(theta_3/2)*exp(-(theta_2*1i)/2)*exp((theta_4*1i)/2);
% M(1,1) = cos(theta_1/2)*cos(theta_3/2)*exp(-(theta_2*1i)/2)*exp(-(theta_4*1i)/2);
theta_4 =  angle(M(2,2)) - angle (M (1,1));

% M(1,4) = -sin(theta_1/2)*sin(theta_3/2)*exp((theta_2*1i)/2)*exp((theta_4*1i)/2);
% M(2,3) = -sin(theta_1/2)*sin(theta_3/2)*exp((theta_2*1i)/2)*exp(-(theta_4*1i)/2);
theta_2 = angle(M(1,4)) + angle(M(2,3));

%Angle correction to [-pi/2, pi/2]
if theta_2 >= 3*pi/4
    theta_2 = 2*pi - theta_2;
end
if theta_2 <= -3*pi/4
    theta_2 = 2*pi + theta_2;
end


% Angle correction to [-pi/2,pi/2]
% M(1,2)= -cos(theta_1/2)*sin(theta_3/2)*exp(-(theta_2*1i)/2)*exp((theta_4*1i)/2)*1i;
if  M(1,2)/(exp(-(theta_2*1i)/2)*exp((theta_4*1i)/2)*1i) >= 0.0
    if theta_3 >= 0.0
       theta_3 = -theta_3;
    end
end

% Angle correction to [-pi/2,pi/2]
% M(1,3) = -cos(theta_3/2)*sin(theta_1/2)*exp((theta_2*1i)/2)*exp(-(theta_4*1i)/2)*1i;
if M(1,3)/( exp((theta_2*1i)/2)*exp(-(theta_4*1i)/2)*1i ) >= 0.0
    if theta_1 >= 0.0
       theta_1 = -theta_1;
    end
end
    
% M(1,1) = cos(theta_1/2)*cos(theta_3/2)*exp(-(theta_2*1i)/2)*exp(-(theta_4*1i)/2); 
if abs(M(1,1) - cos(theta_1/2)*cos(theta_3/2)*exp(-(theta_2*1i)/2)*exp(-(theta_4*1i)/2)) >= 0.001
    theta_2 = -theta_2;
end

%Return results
angles = [theta_1 theta_2 theta_3 theta_4 ];


end



%-------------------------------------------------------------------------
%  D : returns the product of the exponentials of the three square matrices
%      A1, A2, A3 and A4  as  D = exp(-ia1A1/2)*exp(-ia2A2/2)*exp(-ia3A3/2)
%      *exp(-ia4A4/2)
%-------------------------------------------------------------------------
function M = D (angles, A1, A2, A3, A4)
    % Use of matrix exponentials 
    M =  expm(-1i*angles(1)*A1/2)...
        *expm(-1i*angles(2)*A2/2)...
        *expm(-1i*angles(3)*A3/2)...
        *expm(-1i*angles(4)*A4/2);
end
