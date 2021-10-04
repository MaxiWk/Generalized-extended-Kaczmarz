    function medianplot_array = plot_minmax_median_quantiles(mins,maxs,medians,quant25s,quant75s,choose_logy,method_array,iter_save,maxiter,minmaxcolor_dict,quantcolor_dict,linecolor_dict,displayname_dict)

      medianplot_array = zeros(1, length(method_array));  % for legend
      num_iter_array = 1:iter_save:maxiter;
      num_methods = length(method_array);

      for i = 1:num_methods
          minmaxcolor_i = minmaxcolor_dict(method_array{i});

          if choose_logy
            h = fill([num_iter_array  fliplr(num_iter_array)], [log10(maxs(:,i)')  fliplr(log10(mins(:,i))')], minmaxcolor_i,'EdgeColor', 'none');
            set(h,'facealpha', .5)
            quantcolor_i = quantcolor_dict(method_array{i});
            h = fill([num_iter_array  fliplr(num_iter_array)], [log10(quant75s(:,i)')  fliplr(log10(quant25s(:,i))')], quantcolor_i,'EdgeColor', 'none');
            set(h,'facealpha', .5)
            medianplot_array(i) = plot( num_iter_array,log10(medians(:,i)),linecolor_dict(method_array{i}),'LineWidth',2,...
                       'DisplayName',displayname_dict(method_array{i}) );
            ylabel('(log scale)')
          else 
            h = fill([num_iter_array  fliplr(num_iter_array)], [maxs(:,i)'  fliplr(mins(:,i)')], minmaxcolor_i,'EdgeColor', 'none');
            set(h,'facealpha', .5)
            quantcolor_i = quantcolor_dict(method_array{i});
            h = fill([num_iter_array  fliplr(num_iter_array)], [quant75s(:,i)'  fliplr(quant25s(:,i)')], quantcolor_i,'EdgeColor', 'none');
            set(h,'facealpha', .5)
            medianplot_array(i) = plot( num_iter_array,medians(:,i),linecolor_dict(method_array{i}),'LineWidth',2,...
                       'DisplayName',displayname_dict(method_array{i}) );          
          end
          % use 10^ notation for iterations    
          %xt = get(gca, 'xtick');
          %set('xticklabel', sprintf('%1.1e|', xt));

      end
    end