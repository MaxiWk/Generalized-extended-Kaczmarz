function medianplot_array = plot_stopped_quantity(quantity,iterstop_list,choose_logy,method_array,method_counter,iter_save,maxiter,linecolor_dict,displayname_dict)
      
      medianplot_array = zeros(1, length(method_array));  % for legend
      num_iter_array = 1:iter_save:maxiter;
      num_methods = length(method_array);
      
      iter_stop = iterstop_list(end,method_counter);
      index_stop = ceil( iter_stop/iter_save );

      for method_counter = 1:num_methods
          
            if choose_logy
                plotted_quantity = log10(quantity(:,end,method_counter));
                stopped_quantity = log10(quantity(index_stop,end,method_counter));
                ylabel('log scale')
            else
                plotted_quantity = quantity(:,end,method_counter);
                stopped_quantity = quantity(index_stop,end,method_counter);
            end

            medianplot_array(method_counter) = plot(num_iter_array, plotted_quantity', ...
                linecolor_dict(method_array{method_counter}),'LineWidth',2,...
                'DisplayName',displayname_dict(method_array{method_counter}));

            if any( strcmp(method_array{method_counter},{'grek','egrek'}) )
                plot(iter_stop, stopped_quantity, 'o', 'markersize', 10, 'color', 'k');
            end

      end

end

