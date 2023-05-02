function [opt] = telescope_opt(opt,tel_i,output)
%Telescope_opt written by Sarah Palmer Jan 2023
%  This function sets up the grid arrays for each telescoping grid
j = opt.bf.j;
k = opt.bf.k;
l = opt.bf.l;
m = opt.bf.m;
n = opt.bf.n;
pd = opt.pd;
pm = opt.pm;
if pd == 5
    if tel_i == 1
        opt.dies.kW{tel_i} = linspace(opt.dies.kW_1,opt.dies.kW_m,j);              %[kW]
        opt.inso.kW{tel_i} = linspace(opt.inso.kW_1,opt.inso.kW_m,k);              %[kW]
        opt.wind.kW{tel_i} = linspace(opt.wind.kW_1,opt.wind.kW_m,l);              %[kW]
        opt.wave.kW{tel_i} = linspace(opt.wave.kW_1,opt.wave.kW_m,m);              %[kW]
        opt.Smax{tel_i} = linspace(opt.Smax_1,opt.Smax_n,n);                       %[kWh]
    else
        %find the index at the min value from the last telescope
        %opt.dies.kW{tel_i-1} == output.min.kWd{tel_i-1}
        j_min = find(opt.dies.kW{tel_i-1} == output.min.kWd{tel_i-1})
        k_min = find(opt.inso.kW{tel_i-1} == output.min.kWi{tel_i-1})
        l_min = find(opt.wind.kW{tel_i-1} == output.min.kWwi{tel_i-1})
        m_min = find(opt.wave.kW{tel_i-1} == output.min.kWwa{tel_i-1})
        n_min = find(opt.Smax{tel_i-1} == output.min.Smax{tel_i-1})
        
        %find new min and max points of next telescoped grid (accounts for a
        %min at the edge of the tele_i-1 grid
        [j_1,j_max] = check_edges(j, j_min);
        [k_1,k_max] = check_edges(k, k_min);
        [l_1,l_max] = check_edges(l, l_min);
        [m_1,m_max] = check_edges(m, m_min);
        [n_1,n_max] = check_edges(n, n_min);
    
        %Make new arrays
        opt.dies.kW{tel_i} = linspace(opt.dies.kW{tel_i-1}(j_1),opt.dies.kW{tel_i-1}(j_max),j);             %[kW]
        opt.inso.kW{tel_i} = linspace(opt.inso.kW{tel_i-1}(k_1),opt.inso.kW{tel_i-1}(k_max),k);              %[kW]
        opt.wind.kW{tel_i} = linspace(opt.wind.kW{tel_i-1}(l_1),opt.wind.kW{tel_i-1}(l_max),l);              %[kW]
        opt.wave.kW{tel_i} = linspace(opt.wave.kW{tel_i-1}(m_1),opt.wave.kW{tel_i-1}(m_max),m);   
        opt.Smax{tel_i} = linspace(opt.Smax{tel_i-1}(n_1),opt.Smax{tel_i-1}(n_max),n);
    end
elseif pd == 2
    %1:Wi 2:In 3:Wa 4:Di
    if tel_i == 1
        opt.Smax{tel_i} = linspace(opt.Smax_1,opt.Smax_n,n);                       %[kWh]
        disp('test - past battery')
        if pm == 4 
            opt.dies.kW{tel_i} = linspace(opt.dies.kW_1,opt.dies.kW_m,j);              %[kW]
        elseif pm == 2
            opt.inso.kW{tel_i} = linspace(opt.inso.kW_1,opt.inso.kW_m,k);              %[kW]
        elseif pm == 1
            opt.wind.kW{tel_i} = linspace(opt.wind.kW_1,opt.wind.kW_m,l);              %[kW]
        elseif pm == 3
            opt.wave.kW{tel_i} = linspace(opt.wave.kW_1,opt.wave.kW_m,m);              %[kW]
        end
    else
        %find the index at the min value from the last telescope
        %opt.dies.kW{tel_i-1} == output.min.kWd{tel_i-1}
        n_min = find(opt.Smax{tel_i-1} == output.min.Smax{tel_i-1})
        [n_1,n_max] = check_edges(n, n_min);
        opt.Smax{tel_i} = linspace(opt.Smax{tel_i-1}(n_1),opt.Smax{tel_i-1}(n_max),n);
        if pm == 4 
            j_min = find(opt.dies.kW{tel_i-1} == output.min.kWd{tel_i-1})
            [j_1,j_max] = check_edges(j, j_min);
            opt.dies.kW{tel_i} = linspace(opt.dies.kW{tel_i-1}(j_1),opt.dies.kW{tel_i-1}(j_max),j);     
        elseif pm == 2
            k_min = find(opt.inso.kW{tel_i-1} == output.min.kWi{tel_i-1})
            [k_1,k_max] = check_edges(k, k_min);
            opt.inso.kW{tel_i} = linspace(opt.inso.kW{tel_i-1}(k_1),opt.inso.kW{tel_i-1}(k_max),k); 
        elseif pm == 1
            l_min = find(opt.wind.kW{tel_i-1} == output.min.kWwi{tel_i-1})
            [l_1,l_max] = check_edges(l, l_min);
            opt.wind.kW{tel_i} = linspace(opt.wind.kW{tel_i-1}(l_1),opt.wind.kW{tel_i-1}(l_max),l); 
        elseif pm == 3
            m_min = find(opt.wave.kW{tel_i-1} == output.min.kWwa{tel_i-1})
            [m_1,m_max] = check_edges(m, m_min);
            opt.wave.kW{tel_i} = linspace(opt.wave.kW{tel_i-1}(m_1),opt.wave.kW{tel_i-1}(m_max),m);  
        end
        
    end
end
end