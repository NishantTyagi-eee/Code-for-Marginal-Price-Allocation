clear;
clc;

% Define symbolic variables
syms s TR zeta w_n R D Km t t1 t2 beta;
t_values = 0:0.01:35;  % Defining t-axis values
%% STEP 1: Finding the system parameters from the given SFR model 
% Input the SFR transfer function
% Actual_Num=input('Please enter the numerator of SFR\n');
% Actual_Den=input('Please enter the denominator of SFR\n');
% Actual_Const=input('Please enter the constant of SFR\n');
Actual_Num = [3.39 1];
Actual_Den = [1 0.38 0.44^2];
Actual_Const=0.018;
GSFR_Std = ((1+TR*s)/(s^2 + 2*zeta*w_n*s + w_n^2)); % Define the STANDARD G_SFR equation
[Std_Num, Std_Den] = numden(GSFR_Std); % Extract the numerators and denominators of GSFR_Std
Std_Num_Coeffs=coeffs(Std_Num,s); % Extract coefficients of numerator of GSFR_Std
Std_Den_Coeffs=coeffs(Std_Den,s); % Extract coefficients of denominator of GSFR_Std

% Equate coefficients and solve for TR zeta and w_n
% CAUTION: coeffs() function returns coefficients in ascending powers of 's'
% but user entered Num and Den contain coefficients in decreasing powers of 's'
equations = [Std_Num_Coeffs(2) == Actual_Num(1), Std_Num_Coeffs(1) == Actual_Num(2), ...
    Std_Den_Coeffs(3) == Actual_Den(1), Std_Den_Coeffs(2) == Actual_Den(2), Std_Den_Coeffs(1) == Actual_Den(3)];
variables = [TR zeta w_n];
solution = solve(equations, variables);

% Extract ONLY UNIQUE and POSITIVE solution
TR=double(unique(solution.TR));
zeta=double(solution.zeta(solution.zeta > 0));
w_n=double(solution.w_n(solution.w_n > 0));
% "solution.x > 0" creates a logical array that is true for elements of solution.x that are greater than 0. 
% solution.x(solution.x > 0) then uses this logical array to return only the positive solutions.
GSFR_Act_Integral = double(Actual_Const/(w_n^2));

%% STEP 2: Determine the CURRENT frequency response
% Use TR, zeta, w_n to calculate alpha, w_r and phi
alpha = sqrt((1-2*TR*zeta*w_n+(TR*w_n)^2)/(1-zeta^2)); % Calculate alpha
w_r = w_n*sqrt(1-zeta^2); % Calculate w_r
phi= pi + atan((sqrt(1-zeta^2))/(zeta-(TR*w_n))); % Calculate phi
% P = input('Please enter the size of step disturbance in pu\n');
P = 2; % Define the power disturbance to be a step (in pu)
f_0=50; % Nominal frequency
f_t = f_0 - (GSFR_Act_Integral*P) * (1 - alpha * exp(-zeta*w_n*t) * sin(w_r*t + phi)); % Define function f(t)
f_t_values = double(subs(f_t, t, t_values)); % Finding f_t over t-axis
figure(1)
% width = 5; % inches
% height = 3.5; % inches
% figure('Units', 'inches', 'Position', [0, 0, width, height]);
plot(t_values,(f_t_values),'-','LineWidth',2.2,'Color','[0.372, 0.718, 0.631]') % Plotting f_t
hold on

%% STEP 3: Determine the f_f, f_n, and f_n_new:
f_f = f_0 - (GSFR_Act_Integral*P); % Compute the quasi steady state settling frequency, f_f
t_n = (atan((w_r*TR)/(zeta*w_n*TR-1))+pi)/w_r; % Find the nadir time, t_n
f_n = double(subs(f_t, t, t_n)); % Find f_nadir i.e., f(n), by evaluating f(t) at t_n
D_tr = double(f_f - f_n); % Find transient frequency deviation
% Assume a portion beta of transient deviation D_tr is to be removed
% beta=input('Please enter the fraction of transient deviation which is desired to be removed\n');
beta=1;
f_n_req = double(beta*f_f + (1-beta)*f_n); % Find new i.e., desired f_nadir
% yline(f_n_req, 'Color', [0.9, 0.6, 0], 'LineWidth', 1, 'Label', '$f_n^{req}$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12,'FontWeight', 'bold', 'HandleVisibility', 'off');


%% STEP 4: Determine the DESIRED frequency response f_prime(t)
% Finding roots t1 and t2, where curve f(t) intersects the horizontal line given by f_n_new
eqn = f_t == f_n_req;  % Creating the equation to solve
init_guesses = [t_n-1, t_n+1]; % Creating a vector of initial guesses where roots t1 and t2 may lie
t_solutions = zeros(size(init_guesses)); % Declaring a vector to store the roots t1 and t2
for i = 1:length(init_guesses) % Solve the equation for variable t with initial guesses vector
    t_solutions(i) = double(vpasolve(eqn,t,init_guesses(i)));
end
t_solutions = sort(t_solutions(t_solutions > 0)); % t_solutions will contain all roots,
% both positive and negative, we will select only the positive roots and sort them
% t1 and t2 will be the first two positive roots of the equation
t1 = t_solutions(1); % Assign first root to t1
t2 = t_solutions(2); % Assign second root to t2.
f_prime_t = piecewise(t < t1 | t > t2, f_t, t1 <= t & t <= t2, f_n_req); % Define f_prime using t1 and t2
f_prime_t_values = double(subs(f_prime_t, t, t_values)); % Finding f_prime over t-axis
figure(1)
plot(t_values,(f_prime_t_values),':','LineWidth',2.2,'Color','[0.543, 0.169, 0.886]') % Plotting f_prime
hold on
%% STEP 5: Determine the CHANGE in frequency response, ∆f(t):
delta_f_t = (f_prime_t - f_t); 
delta_f_t_values = double(subs(delta_f_t, t, t_values)); % Evaluate Δf(t) over t-axis
plot(t_values,(delta_f_t_values + 50),'-','LineWidth',2.2,'Color','[0.802, 0.343, 0.294]');
lgd = legend('$f^{act}(t)$', '$f^{req}(t)$ with $\beta=1$', '$\Delta f^{comp}(t)$');
set(lgd, 'Interpreter', 'latex', 'FontName', 'Times New Roman','Fontsize',12);
xlabel('Time (s)', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12);
set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontName', 'Times New Roman');
set(gca, 'YTickLabel', get(gca, 'YTickLabel'), 'FontName', 'Times New Roman');

hold on;
grid on;

%% STEP 6: Extract the last i.e., reqiured section of Δf(t)'s expression
delta_f_t_str = char(delta_f_t); % Convert symbolic piecewise function to string
parts = strsplit(delta_f_t_str, ',');% Split string into parts by comma
required_part_str = strtrim(parts{end}(1:end-1));% Required part of Δf(t)'s expression is the second last element in "parts", excluding the last character ")"
delta_f_t_sym = str2sym(required_part_str); % We need Δf(t)'s expression as a symbol instead of a string
% Define rectangular pulse between t1 and t2 so that Δf(t) is considered only in [t1, t2]
syms pulse
pulse = heaviside(t - t1) - heaviside(t - t2);
delta_f_t_sym=delta_f_t_sym*pulse;

%% STEP 7: Determine Δp(t)
GSFR_Actual = Actual_Const * (poly2sym(Actual_Num, s) / poly2sym(Actual_Den, s)); % Construct GSFR_Actual
delta_f_s=laplace(delta_f_t_sym); % Find Laplace transform of Δf(t)
delta_p_s = delta_f_s/GSFR_Actual; % Now, find ΔPs = ΔFs/G_SFR(s)
delta_p_t=ilaplace(delta_p_s); % Last, find Δp(t)
delta_p_t_values = double(subs(delta_p_t, t, t_values));
figure(2)
% width = 5; % inches
% height = 3.5; % inches
% figure('Units', 'inches', 'Position', [0, 0, width, height]);
plot(t_values,(delta_p_t_values*800),'-','LineWidth',2.2,'Color','[0.802, 0.343, 0.294]');
hold on
grid on

%% STEP 8: Finding the equivalent triangular Δp(t) with (a) equal area, (b) equal height, (c) same start time, (d) rise time of t_rise seconds
tri_start = t1-0.2; % start instant of original Δp(t) and equivalent triangular Δp(t) is same
max_height = max(delta_p_t_values); % Calculate max height of the original curve
time_end=25; % Calculate the time range for integrating the original curve
% end time is a reasonable upper limit where we know the curve has decayed
delta_p_t_area = double(int(delta_p_t, t1, time_end)); % calculate area under the original curve
% for a triangle, area = 0.5*base*height ==> calculate base of the equivalent triangle
tri_base = (2 * delta_p_t_area) / max_height;
tri_end=tri_start+tri_base;
% Define equivalent triangular Δp(t) function using symbolic piecewise function
tri_rise=0.5;
tri_peak=tri_start+tri_rise;
tri_fall=(tri_base+tri_start)-tri_peak;
% slope_rising = max_height / tri_rise; % Calculate slope of rising edge
% slope_falling = -max_height / tri_fall; % Calculate slope of falling edge
rising_edge = (max_height/tri_rise) * (t - tri_start);
falling_edge =  - ((max_height/tri_fall) * (t - tri_end));
tri_delta_p_t = piecewise(t < tri_start, 0, tri_start <= t & t < tri_peak, rising_edge, tri_peak <= t & t <= tri_end, falling_edge, t > tri_end, 0);
% Plot the triangular Δp(t)
tri_delta_p_t_values=double(subs(tri_delta_p_t,t,t_values));
figure(2)
plot(t_values,tri_delta_p_t_values*800,'-.','LineWidth',2.2,'Color','[0.274, 0.509, 0.705]');
lgd = legend('$\Delta p^{comp}(t)$', '$ \Delta p^{comp}_{tri}(t)$');
set(lgd, 'Interpreter', 'latex', 'FontName', 'Times New Roman','Fontsize',12);
xlabel('Time (s)', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Power Injection (MW)', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12);
set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontName', 'Times New Roman');
set(gca, 'YTickLabel', get(gca, 'YTickLabel'), 'FontName', 'Times New Roman');


%% STEP 9: Finding the new Δf(t) when new equivalent triangular Δp(t)
tri_delta_p_t_rise = rising_edge*(heaviside(t-tri_start)-heaviside(t-tri_peak));
tri_delta_p_t_fall=falling_edge*(heaviside(t-tri_peak)-heaviside(t-tri_end));
tri_delta_p_s_rise=laplace(tri_delta_p_t_rise);
tri_delta_p_s_fall=laplace(tri_delta_p_t_fall);
tri_delta_p_s=tri_delta_p_s_rise+tri_delta_p_s_fall;

t_start_delta_p_t_2=ceil(tri_end);
t_end_delta_p_t_2=t_start_delta_p_t_2+10;
delta_p_t_2_rise=0.1*(t-10)*(heaviside(t-9)-heaviside(t-20));
delta_p_t_2_constant=1*heaviside(t-20);
delta_p_s_2_rise=laplace(delta_p_t_2_rise);
delta_p_s_2_constant=laplace(delta_p_t_2_constant);
delta_p_s_2=delta_p_s_2_rise+delta_p_s_2_constant;

%delta_p_s_total = tri_delta_p_s + delta_p_s_2;

% delta_f_s_new_1=tri_delta_p_s*GSFR_Actual;
% delta_f_t_new_1=ilaplace(delta_f_s_new_1);

% delta_f_s_new_2=delta_p_s_2*GSFR_Actual;
% delta_f_t_new_2=ilaplace(delta_f_s_new_2);
delta_p_s_new=tri_delta_p_s+delta_p_s_2;
delta_f_s_new=delta_p_s_new*GSFR_Actual;
delta_f_t_new=ilaplace(delta_f_s_new);
f_t_new=f_t+delta_f_t_new;
f_t_new_values=double(subs(f_t_new,t,t_values));
figure(3)
ax = gca;
ax.YLim = [49.7 50.1];
yline(50, 'Color', 'k', 'LineWidth', 1, 'Label', '$f_0$', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12,'FontWeight', 'bold', 'HandleVisibility', 'off');
hold on
plot(t_values,(f_t_values),'-','LineWidth',2.2,'Color','[0.372, 0.718, 0.631]') % Plotting f_t
hold on
plot(t_values,(f_t_new_values),'-.','LineWidth',2.2,'Color','[0.543, 0.169, 0.886]')
hold on
lgd = legend('$f^{act}(t)$ after LoG', '$f^{res}(t)$ after proposed power injection');
set(lgd, 'Interpreter', 'latex', 'FontName', 'Times New Roman','Fontsize',11);
xlabel('Time (s)', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12);
set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontName', 'Times New Roman');
set(gca, 'YTickLabel', get(gca, 'YTickLabel'), 'FontName', 'Times New Roman');
box on



figure(4)
FFR1_tri_delta_p_t = piecewise(tri_start <= t & t < tri_peak, rising_edge, tri_peak <= t & t <= tri_end, falling_edge);
% Plot the triangular Δp(t)
FFR1_tri_delta_p_t_values=double(subs(FFR1_tri_delta_p_t,t,t_values));
plot(t_values,FFR1_tri_delta_p_t_values*800,'-.','LineWidth',2.2,'Color','[0.274, 0.509, 0.705]');
hold on
FFR2_delta_p_t = piecewise(10 <= t & t < 20, delta_p_t_2_rise, t > 20, delta_p_t_2_constant);
FFR2_delta_p_t_values=double(subs(FFR2_delta_p_t,t_values));
plot(t_values,(FFR2_delta_p_t_values*800),'-','LineWidth',2.2,'Color','[0.543, 0.169, 0.886]')
hold on
lgd = legend('$FR_{1}(t)$ addresses Transient Deviation', '$FR_{2}(t)$ addresses QSS Deviation'); 
set(lgd, 'Interpreter', 'latex', 'FontName', 'Times New Roman','Fontsize',11);
xlabel('Time (s)', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12);
ylabel('Proposed Power Injection (MW)', 'Interpreter', 'latex', 'FontName', 'Times New Roman', 'FontSize', 12);
set(gca, 'XTickLabel', get(gca, 'XTickLabel'), 'FontName', 'Times New Roman');
set(gca, 'YTickLabel', get(gca, 'YTickLabel'), 'FontName', 'Times New Roman');


hold on
