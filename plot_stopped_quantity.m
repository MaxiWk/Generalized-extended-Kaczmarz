function plot_array = plot_stopped_quantity(quantity,iterstop_list,choose_logy,method_array,iter_save,maxiter,linecolor_dict,displayname_dict)

plot_array = zeros(1, length(method_array));  % for legend

for i = 1:length(method_array)
        
        iter_stop = iterstop_list(end,i);
        index_stop = ceil( iter_stop/iter_save ); 
        num_iter_array = 1:iter_save:maxiter;

        if choose_logy
            plot_array(i) = plot(iter_stop, log10(quantity(index_stop,end,i)), '.', 'markersize',30);
            hold on
            plot_array(i) = plot(num_iter_array, log10(quantity(:,end,i)'), ...
                'LineWidth',2,'DisplayName',displayname_dict(method_array{i}));
                        
        else
            plot_array(i) = plot(iter_stop, quantity(index_stop,end,i), '.', 'markersize',30);
            hold on
            plot_array(i) = plot(num_iter_array, quantity(:,end,i)', ...
                'LineWidth',2,'DisplayName',displayname_dict(method_array{i}));                        
        end            
        
end


end