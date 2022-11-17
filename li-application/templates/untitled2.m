clear
rng(1);
f = fopen('T1_wgs_state_clocks.xml');
g = fopen('T1_wgs_state_clocks_red.xml','w');
first = true;
while ~feof(f)
    line = fgets(f);
    if contains(line, '<sequence')
        tmp = strsplit(line, '"');
        if first
            ind = sort(randsample(length(tmp{10}), 1000,1));
            first = false;
        end
        fprintf(g, strrep(line, tmp{10}, tmp{10}(ind)));
    else
        fprintf(g, line)
    end
end
fclose('all');

clear
rng(1);
f = fopen('T2_wgs_state_clocks.xml');
g = fopen('T2_wgs_state_clocks_red.xml','w');
first = true;
while ~feof(f)
    line = fgets(f);
    if contains(line, '<sequence')
        tmp = strsplit(line, '"');
        if first
            ind = sort(randsample(length(tmp{10}), 1000,1));
            first = false;
        end
        fprintf(g, strrep(line, tmp{10}, tmp{10}(ind)));
    else
        fprintf(g, line)
    end
end
fclose('all');




f = fopen('T1_state_clocks_wgs_amplicon_300_sites_2.xml');
g = fopen('T1_state_clocks_wgs_amplicon_300_sites_2_red.xml','w');
first = true;
while ~feof(f)
    line = fgets(f);
    if contains(line, '<sequence')
        tmp = strsplit(line, '"');
        if first
            ind = find(tmp{10}~='N');
            ind = sort([ind, randsample(find(tmp{10}=='N'), 2000,1)]);
            first = false;
        end
        fprintf(g, strrep(line, tmp{10}, tmp{10}(ind)));
    else
        fprintf(g, line)
    end
end
fclose('all');

