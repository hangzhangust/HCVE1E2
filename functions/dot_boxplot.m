function [x_data,y_data] = dot_boxplot(data,nbins,center,max_range,scale,num_interval)
% 
% data = all_E(setdiff(1:L,polymorphisms_associated_with_neutralization_resistance - 383));
% nbins = 50;
% 
% 
% max_range = 0.2;
% 
% 
% center = 1;

[c,b] = hist(data,nbins);
x_data=[];
y_data=[];
data_interval = 2*max_range/num_interval;

c = floor(c/scale);
I = c>num_interval;
c(I)=num_interval;

for i =1:length(b)
    num_points = c(i);
    point_side = floor(num_points/2);
    
    
    
    if num_points>1
%          if  mod(num_points,2)==0
            for j =1: point_side
                y_data = [y_data; b(i)];

                x_data = [x_data;center-point_side*data_interval+(j-1)*data_interval+data_interval/2];
            end

            for j =point_side+1:num_points
                y_data = [y_data; b(i)];

                x_data = [x_data;center-point_side*data_interval+j*data_interval-data_interval/2];
            end
%          end
    else
        
        if  mod(num_points,2)==1
            y_data = [y_data; b(i)];
            x_data = [x_data; center];
             
        end
        
    end
end
end