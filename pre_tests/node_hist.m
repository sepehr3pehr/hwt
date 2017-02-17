nb = 128;

hist_fl = zeros(nb,1);
lookup = zeros(256,1);


leve = 1;
num_data = 10^9;
%data = load();

nb_chunk = nb/8;

for i=1:256
    t = dec2bin(i-1);
    lookup(i) = length(find(t=='1'));
end


for i=1:num_data
    if(mod(i,10000)==0)
        fprintf('i = %d',i);
    end
    ind = sum(lookup(data(:,i)));
    ind = ind + 1;
    hist_fl(ind) = hist_fl(ind) + 1;
end 
