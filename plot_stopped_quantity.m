function plot_stopped_quantity(quantity,iterstop_list,choose_logy,method_array,method_counter,iter_save,maxiter)
      
    iter_stop = iterstop_list(end,method_counter);
    index_stop = ceil( iter_stop/iter_save ); 
    num_iter_array = 1:iter_save:maxiter;

    if choose_logy
        plot(iter_stop, log10(quantity(index_stop,end,method_counter)), '.', 'markersize', 30, 'color', 'red');
        hold on
        plot(num_iter_array, log10(quantity(:,end,method_counter)'), ...
            'LineWidth',2,'DisplayName',method_array{method_counter}, 'color', 'blue');

    else
        plot(iter_stop, quantity(index_stop,end,method_counter), '.', 'markersize', 30, 'color', 'red');
        hold on        
        plot(num_iter_array, quantity(:,end,method_counter)', ...
            'LineWidth',2,'DisplayName',method_array{method_counter}, 'color', 'blue');                        
    end            

    

end

