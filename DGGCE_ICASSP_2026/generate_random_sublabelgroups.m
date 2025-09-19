function [sublabelgroups, BPi_sublabelgroups] = generate_random_sublabelgroups(BPi, K, m, seed)

    rng(seed);  
    
    num_labels = size(BPi, 2);  
    all_indices = 1:num_labels;  
    
    sublabelgroups = cell(K, 1);  
    BPi_sublabelgroups = cell(K, 1); 
    used_labels = false(1, num_labels);  
    

    for i = 1:num_labels
        chosen_subspace = randi(K);  
        sublabelgroups{chosen_subspace} = [sublabelgroups{chosen_subspace}, i]; 
    end
    

    for k = 1:K

        while length(sublabelgroups{k}) < m
           
            random_label = randi(num_labels);  
            sublabelgroups{k} = [sublabelgroups{k}, random_label]; 
        end
   
        BPi_sublabelgroups{k} = BPi(:, sublabelgroups{k});
    end
end