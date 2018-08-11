import heapq


def generate_ss_configuration_freq(most_frequent_1,total_dict_1,complete_total_1,ss_1,most_frequent_2,total_dict_2,complete_total_2,ss_2):
    ss_title_1 = str(ss_1).replace(')','')
    ss_title_1 = ss_title_1.replace('(','')
    ss_title_2 = str(ss_2).replace(')','')
    ss_title_2 = ss_title_2.replace('(','')
    fig = plt.subplot(121)
    #fig, ax = plt.subplots()
    i= 0
    config_ticks = []
    for key in most_frequent_1:
        config = total_dict_1[key]
        config = str(config).replace(')','')
        config = config.replace('(','')
        config_ticks.append(config)
        try:
            plt.bar(i, float(key)/float(complete_total_1),color = 'blue')
        except ZeroDivisionError:
            plt.bar(i, 0,color = 'blue')
        i =i+1
    plt.ylim(0,1.23)
    xticks_pos = [i for i in range(0,len(config_ticks))]
    plt.xticks(xticks_pos,x_axis_labels,rotation=45)
    plt.title(ss_title_1)
    plt.ylabel("Frequency")
    plt.xlabel('Configuration')

    fig = plt.subplot(122)
    #fig, ax = plt.subplots()
    i= 0
    config_ticks = []
    for key in most_frequent_2:
        config = total_dict_2[key]
        config = str(config).replace(')','')
        config = config.replace('(','')
        config_ticks.append(config)
        try:
            plt.bar(i, float(key)/float(complete_total_2),color = 'blue')
        except ZeroDivisionError:
            plt.bar(i, 0,color = 'blue')
        i =i+1
    plt.ylim(0,1.23)
    xticks_pos = [i for i in range(0,len(config_ticks))]
    plt.xticks(xticks_pos,x_axis_labels,rotation=45)
    plt.title(ss_title_2)
    plt.ylabel("Frequency")
    plt.xlabel('Configuration')


    plt.show()
    return()


df['combined'] = df[['Cys1_SS_cat','Cys2_SS_cat']].apply(lambda x: ','.join(x), axis=1)

k = 0
for ss in unique_ss:
    ss = ss.split(',')
    ss_forward = ss[0] + ',' +ss[1]
    ss_reverse = ss[1]+ ',' +ss[0]

    if ss_forward == ss_reverse:
        ss_df = df.loc[df['combined'] == ss_forward]        
    if ss_forward != ss_reverse:
        ss_df = df.loc[(df['combined'] == ss_forward) | (df['combined'] == ss_reverse)]

    total = len(ss_df)
    if total > 100:
        if k == 0:

            ss_1 = ss
            total_dict_1 = {}
            total_list_1 = []
            complete_total_1 = 0
            
            for config in configurations:
                    common_list=[]
                    config_total = 0
                    forward_config = 0
                    reverse_config = 0
                    forward_config = len(ss_df.loc[
                            (ss_df['Cys1_x1'] == float(config[0])) & 
                            (ss_df['Cys1_x2'] == float(config[1])) & 
                            (ss_df['x3'     ] == float(config[2])) & 
                            (ss_df['Cys2_x2'] == float(config[3])) & 
                            (ss_df['Cys2_x1'] == float(config[4])) 
                           ])
        
        ###################################################################################
        # If the configuration IS NOT symmetrical, then search for the reverse configuration
        # x1b,x2b,x3,x2,x1
        ###################################################################################
                    if config != config[::-1]:
                        reverse_config = len(ss_df.loc[
                                 (ss_df['Cys2_x1'] == float(config[0])) & 
                                 (ss_df['Cys2_x2'] == float(config[1])) & 
                                 (ss_df['x3'     ] == float(config[2])) & 
                                 (ss_df['Cys1_x2'] == float(config[3])) & 
                                 (ss_df['Cys1_x1'] == float(config[4])) 
                                ])
                    config_total = forward_config + reverse_config
                    #if config_total > 0:
                    total_list_1.append(config_total)
                    total_dict_1[config_total] = config
                    common_list.append(config_total)
                    
                    complete_total_1 = complete_total_1 + config_total
            most_frequent_1 = heapq.nlargest(5,total_list_1)
       

        if k == 1:
            ss_2 = ss

            total_dict_2 = {}
            total_list_2 = []
            complete_total_2 = 0
            
            for config in configurations:
                    common_list=[]
                    config_total = 0
                    forward_config = 0
                    reverse_config = 0
                    forward_config = len(ss_df.loc[
                            (ss_df['Cys1_x1'] == float(config[0])) & 
                            (ss_df['Cys1_x2'] == float(config[1])) & 
                            (ss_df['x3'     ] == float(config[2])) & 
                            (ss_df['Cys2_x2'] == float(config[3])) & 
                            (ss_df['Cys2_x1'] == float(config[4])) 
                           ])
        
        ###################################################################################
        # If the configuration IS NOT symmetrical, then search for the reverse configuration
        # x1b,x2b,x3,x2,x1
        ###################################################################################
                    if config != config[::-1]:
                        reverse_config = len(ss_df.loc[
                                 (ss_df['Cys2_x1'] == float(config[0])) & 
                                 (ss_df['Cys2_x2'] == float(config[1])) & 
                                 (ss_df['x3'     ] == float(config[2])) & 
                                 (ss_df['Cys1_x2'] == float(config[3])) & 
                                 (ss_df['Cys1_x1'] == float(config[4])) 
                                ])
                    config_total = forward_config + reverse_config
                    #if config_total > 0:
                    total_list_2.append(config_total)
                    total_dict_2[config_total] = config
                    common_list.append(config_total)
                    
                    complete_total_2 = complete_total_2 + config_total
            
        k = 0   
        most_frequent_2 = heapq.nlargest(5,total_list_2
        generate_ss_configuration_freq(most_frequent_1,total_dict_1,complete_total_1,ss_1,most_frequent_2,total_dict_2,complete_total_2,ss_2)
        k = k+1