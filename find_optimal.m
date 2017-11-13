%Parameters: mass (m, kg), acceleration due to gravity (g, m(s^-2)), maximum power
%(p_max, W), max_speed (v_max, ms^-1), max_electric slope (max_s, W), brake
%conversion (c_rate)
%(GOING TO ADD LATER) maximum charge (c_max, J or MWH), minimum charge
%(c_min J or MWH), regeneration rate (r, % of power), fuel conversion parameters (d_conv, e_conv, %of power), parameters (c_min J or MWH)

m=13409; %kg, hybrid; regular 11364

m_regular=11364;

g=9.8; %ms^-2

p_max=290*10^3;%Watts, max power that can be outputted total

a_max=0; 

v_max=48000/3600;%meters per second

params=[m,g,p_max,a_max,v_max];
%engine parameters

c_max=25*10^3; %kwh
c_min=0.6*c_max; %kwh
d_conv=0.07; 
e_conv=0.23;
max_s=0.0856;
r=0.2;
c_rate=0; 
engine_params=[c_max, c_min, d_conv, e_conv, max_s, r, c_rate];
%Equations: velocity, acceleration, power, position, slope. Position and
%slope are known

%Intiial parameters velocity, accelleration, power, slope and distance of interval at an index (t)
%velocity is initialiazed by velocity after one second of acceleration. 

t=1;
v_t=v_max;
p_t=v_t*m*g*slope(t,1);
a_t=0;


%Final computations for fuel consumption, total fuel used for hybrids are in
%dieselstats(1,:). total fuel used for diesel busses are in
%dieselstats(2,:). Columns represent routes in ascending order. 

dieselstats=zeros(2,6);



[avg_diesel_hybrid, avg_electricity, avg_diesel_diesel]=loop(params,  engine_params, 20, '10route.csv');
dieselstats(1,1)=avg_diesel_hybrid;
dieselstats(2,1)=avg_diesel_diesel;
[h11, e11, d11]=loop(params, engine_params, 20,'11route_day_regular.csv' );
dieselstats(1,2)=h11;
dieselstats(2,2)=d11;
[h15, e15, d15]=loop(params, engine_params, 20,'15route.csv' );
dieselstats(1,3)=h15;
dieselstats(2,3)=d15;    
[h17, e17, d17]=loop(params, engine_params, 20,'17route.csv' );
dieselstats(1,4)=h17;
dieselstats(2,4)=d17;    
[h81, e81, d81]=loop(params, engine_params, 20,'81route.csv' );
dieselstats(1,5)=h81;
dieselstats(2,5)=d81;   
[h82, e82, d82]=loop(params, engine_params, 20,'82route.csv' );
dieselstats(1,6)=h82;
dieselstats(2,6)=d82;    

dieselstats=real(dieselstats)/1000


%Sensitivity/Robustness calculations
sensitivity_matrix=zeros(6,4,2);
testing=zeros(6,5);

testing(:,1)=linspace(0.7*engine_params(1,5),1.3*engine_params(1,5),6).';%max_slope
testing(:,2)=linspace(0.7*engine_params(1,6),1.3*engine_params(1,6),6).';%max_brakeregeneration
testing(:,3)=linspace(0.7*engine_params(1,4),1.3*engine_params(1,4),6).';%battery_efficiency
testing(:,4)=linspace(0.7*engine_params(1,3),1.3*engine_params(1,3),6).';%engine_efficiency
testing(:,5)=[0, 0.001, 0.01, 0.04, 0.08, 0.12];

slope_test=repmat(engine_params,6);
slope_test(:,5)=testing(:,1);
max_regeneration=repmat(engine_params,6);
max_regeneration(:,6)=testing(:,2);
battery_e=repmat(engine_params,6);
battery_e(:,4)=testing(:,3);
engine_e=repmat(engine_params,6);
engine_e(:,3)=testing(:,4);
charging_r=repmat(engine_params,6);
charging_r(:,7)=testing(:,5);

for i=1:6
    [dh, ae, dd]=loop(params,  slope_test(i,:), 20, '10route.csv');
    sensitivity_matrix(i, 1, 1)=dh/dd;
    [dh1, ae1, dd1]=loop(params,  slope_test(i,:), 20, '17route.csv');
    sensitivity_matrix(i, 1, 2)=dh1/dd1;
    [dh2, ae2, dd2]=loop(params,max_regeneration(i,:)  , 20, '10route.csv');
    sensitivity_matrix(i, 2, 1)=dh2/dd2;
    [dh21, ae21, dd21]=loop(params,  max_regeneration(i,:), 20, '17route.csv');
    sensitivity_matrix(i, 2, 2)=dh21/dd21;
    [dh3, ae3, dd3]=loop(params,  battery_e(i,:), 20, '10route.csv');
    sensitivity_matrix(i, 3, 1)=dh3/dd3;
    [dh13, ae13, dd13]=loop(params,  battery_e(i,:), 20, '17route.csv');
    sensitivity_matrix(i, 3, 2)=dh13/dd13;
    [dh4, ae4, dd4]=loop(params,  engine_e(i,:), 20, '10route.csv');
    sensitivity_matrix(i, 4, 1)=dh4/dd4;
    [dh14, ae14, dd14]=loop(params,  engine_e(i,:), 20, '17route.csv');
    sensitivity_matrix(i, 4, 2)=dh14/dd14;
    [dh5, ae5, dd5]=loop(params,  charging_r(i,:), 20, '10route.csv');
    sensitivity_matrix(i, 5, 1)=dh5/dd5;
    [dh15, ae15, dd15]=loop(params,  charging_r(i,:), 20, '17route.csv');
    sensitivity_matrix(i, 5, 2)=dh15/dd15;

end
 
sensitivity_matrix
%function that continues running the bus on a loop


function [avg_diesel_hybrid, avg_electricity, avg_diesel_diesel]=loop(params,  engine_params, num_trials, filename)
%engine parameters
c_max=engine_params(1,1); %kwh
c_min=engine_params(1,2); %kwh
d_conv=engine_params(1,3); 
e_conv=engine_params(1,4);
max_s=engine_params(1,5);%maximum slope
r=engine_params(1,6);
c_rate=engine_params(1,7);

%other parameters
m=params(1,1);

g=params(1,2);

p_max=params(1,3);

a_max=params(1,4);

v_max=params(1,5);


[slope, position, distance]=finds_slope(filename);
slope1=slope;
position1=position;
distance1=distance;
% figure
% plot(slope)
% legend('slope')
t=1;
v_t=v_max;
p_t=v_t*m*g*slope(t,1);
a_t=0;


for trial=1:num_trials
    slope1=vertcat(slope1,slope);
    position1=vertcat(position1,position+position(1,1));
    distance1=vertcat(distance1, distance);
end


[v a p]=progression(v_t, a_t, p_t, slope1, distance1, params);


[diesel, electricity, current]=hybrid_conversion(v, a, p, distance1, slope1, engine_params);
[diesel_mm]=diesel_conversion(v, a, p, distance1, engine_params);

% figure
% plot(diesel)
% legend('diesel')
% figure
% plot(electricity)
% legend('elec')
% figure
% plot(current)
% legend('batt')
avg_diesel_hybrid=diesel(1,end-1)/(num_trials+1);
avg_electricity=electricity(1,end-1)/(num_trials+1);
avg_diesel_diesel=diesel_mm(1,end-1)/num_trials;
end

%function that converts data into fuel used for a diesel truck
function [diesel_mm]=diesel_conversion(v, a, p, distance, engine_params)
%calculate time taken per each step, given distance velocity and
%acceleration
time=zeros(1,length(a));
for i=1:length(a)-1
    time(1,i)=distance(i,1)/((v(1,i)+v(1,i+1))/2);
%     if a(1,i)>0
%        time(1,i)=(-v(1,i)+(v(1,i).^2-4*(1/2)*(a(1,i)).*distance(i,1)).^1/2)./a(1,i);
%     else
%        time(1,i)=distance(i,1)./v(1,i);
%     end
end

d_conv=engine_params(1,3);
c_rate=engine_params(1,7); 


%total diesel used so far
diesel_mm=zeros(1,length(p));%total diesel at each time step

%for each discrete segment
for i=1:length(time)-1
    %if the slope is low enough use electric power
    if p(1,i)>0
        diesel_mm(1,i+1)=diesel_mm(1,i)+(1/d_conv)*p(1,i)*(time(1,i)/3600);
    end
    if p(1,i)<0
        diesel_mm(1,i+1)=diesel_mm(1,i);
    end
end


end

%function that converts power into fuel used, and electric power used for
%hybrid
function [diesel_mm, electricity_mm, c_current_mm, time]=hybrid_conversion(v, a, p, distance, slope, engine_params)
%calculate time taken per each step, given distance velocity and
%acceleration
time=zeros(1,length(a));
for i=1:length(a)-1
    time(1,i)=distance(i,1)/((v(1,i)+v(1,i+1))/2);
%     if a(1,i)>10^-1
%        time(1,i)=(-v(1,i)+(v(1,i).^2-4*(1/2)*(a(1,i)).*distance(i,1)).^1/2)./a(1,i);
%     else
%        time(1,i)=distance(i,1)./v(1,i);
%     end
end
%restating engine parameters

c_max=engine_params(1,1); %kwh
c_min=engine_params(1,2); %kwh
d_conv=engine_params(1,3); 
e_conv=engine_params(1,4);
max_s=engine_params(1,5);%maximum slope
r=engine_params(1,6);
c_rate=engine_params(1,7);

slopes=slope(1:length(slope)-1);

slopes_logical=slopes<max_s;

%total diesel used so far
diesel_mm=zeros(1,length(p));
%total electricity used so far
electricity_mm=zeros(1,length(p));
%total battery charge left
c_current_mm=zeros(1,length(p));
c_current_mm(1,1)=c_max;

%for each discrete segment
for i=1:length(slopes_logical)-1
    %don't allow battery power to go over max power 
    for j=1:length(c_current_mm)-1
        if c_current_mm(1,j+1)>c_max
            c_current_mm(1,j+1)=c_max;
        end
    end
    %if the slope is low enough use electric power
    if slopes_logical(i,1)==1 && c_current_mm(1,i)>c_min
        if p(1,i)>0
            %increment total electricity used, total battery remaining, and
            %total diesel used
            electricity_mm(1,i+1)=electricity_mm(1,i)+(1/e_conv)*p(1,i)*(time(1,i)/3600);
            c_current_mm(1,i+1)=c_current_mm(1,i)-(1/e_conv)*p(1,i)*(time(1,i)/3600);
            diesel_mm(1,i+1)=diesel_mm(1,i);
        end
        %if power outputted is negative charge the battery with brakes 
        if p(1,i)<0
            %increment total electricity used, total battery remaining, and
            %total diesel used
            electricity_mm(1,i+1)=electricity_mm(1,i)+r*p(1,i)*(time(1,i)/3600);
            c_current_mm(1,i+1)=c_current_mm(1,i)-r*p(1,i)*(time(1,i)/3600);
            diesel_mm(1,i+1)=diesel_mm(1,i);
        end
    end
    %if no more battery or slopes are steep 
    if slopes_logical(i,1)==0 || c_current_mm(1,i)<c_min
        %use diesel and charge batterry with engine if power needed is
        %positive
        if p(1,i)>0
            diesel_mm(1,i+1)=diesel_mm(1,i)+(1+c_rate)*(1/d_conv)*p(1,i)*(time(1,i)/3600);
            c_current_mm(1,i+1)=c_current_mm(1,i)+c_rate*(1/d_conv)*p(1,i)*(time(1,i)/3600);
            electricity_mm(1,i+1)=electricity_mm(1,i);
        end
        %charge batterry with engine if power needed is
        %negative
        if p(1,i)<0        
            diesel_mm(1,i+1)=diesel_mm(1,i);
            c_current_mm(1,i+1)=c_current_mm(1,i)-r*p(1,i)*(time(1,i)/3600);
            electricity_mm(1,i+1)=electricity_mm(1,i)+r*p(1,i)*(time(1,i)/3600);
        end
    end
end

end


%function that generates velocity, accelaration, and power information at
%each node of interpolation

function [v a p]=progression(v_t, a_t, p_t, slope, distance, params)
%parameters

m=params(1,1);

g=params(1,2);

p_max=params(1,3);

a_max=params(1,4);

v_max=params(1,5);

%initialize vectors 

v=zeros(1,length(slope));
a=zeros(1,length(slope));
p=zeros(1,length(slope));
v(1,1)=v_max;
p(1,1)=v(1,1)*m*g*slope(1,1);
a(1,1)=0;

for t=1:length(slope)-1
    values=update(v(1,t), a(1,t), p(1,t), t, slope, distance, params);
    v(1,t+1)=values(1,1);
    a(1,t+1)=values(1,2);
    p(1,t+1)=values(1,3);
end


end

%update rule function 


function updated=update(v, a, p, t, slope, distance, params)

%parameters


m=params(1,1);

g=params(1,2);

p_max=params(1,3);

a_max=params(1,4);

v_max=params(1,5);

p1=slope(t,1)*v_max*m*g;
%if power outputted is max because of the slope then maintain max power and decrease velocity 

% if v_max*(slope(t+1,1)*m*g)>=p_max %this will have to change when we add friction
%     p_1=p_max;
%     v_1=p_max/(slope(t+1,1)*m*g); %this will have to change when we add friction
%     a_1=0;
% end

%if power needed to get up the slope is not max then  

% if v_max*(slope(t+1,1)*m*g)<p_max
%     %if velocity is 30 then allow power to vary in order to maintain velocity
%     if (v-v_max)>-1
%         v_1=v_max;
%         a_1=0;
%         p_1=m*g*slope(t+1,1)*v_1;
%     end
%     %if velocity is less than 30 the accelerate using the power not used to
%     %fight gravity
%     % in order to increase velocity to v_max
%     
%     if v-v_max<-1
%         v_1=v+a*(distance(t,1)/v);
%         a_1=((p_max-v*(slope(t,1)*m*g))/((distance(t,1)/v)*2*m))^(1/2);
%         p_1=m*g*slope(t+1,1)*v_1;
%         
%     end
% end

%returns velocity, acceleration and power for that time interval, and next
%time index
updated=[v_max, 0, p1, t+1];


end





function [slope, position, distance] = finds_slope(filename)
    % find slope of data in `filename`
    
    % assumes first row contains header
    M = csvread(filename, 1);
    
    % first column has distances and second column has corresponding
    % elevation
    
    M2 = M(1:end - 1, :);
    M3 = M(2:end, :);
    
    delM = M3 - M2;
    
    % del elevation / del distance
    %slope = smoothdata(delM(:, 2), 'gaussian', 20) ./ smoothdata(delM(:, 1),'gaussian', 20);
    slope=delM(:,2)./delM(:,1);
    position=M(:,1);
    distance=(delM(:,1).^2+delM(:,2).^2).^1/2;
    % comment out if necessary
    %plot(M(:, 1), M(:, 2))
end