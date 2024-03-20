def dec2bin(number,n_site):
    s = [0]*n_site
    index=0
    while number > 0:
        rem = number % 2
        if rem==1:
            s[index]=rem
        index+=1
        number = number // 2
    return s

def bin2dec(bin_list):
    dec=0
    for i in range(len(bin_list)):
        if bin_list[i]==1:
            dec+=2**i
    return dec

def get_reversed_dec(dec,n_site):
    return bin2dec(list(reversed(dec2bin(dec,n_site))))

def get_inverse_amp_indices(n_qubit):
    dim=2**n_qubit
    mapper=[0]*dim
    for amp_index in range(dim):
        mapper[amp_index]=get_reversed_dec(amp_index,n_qubit)
    return mapper

def get_mapped_amp_vec(n_qubit,amp_vec):
    dim=n_qubit**2
    mapper=get_inverse_amp_indices(n_qubit)
    new_amp_vec=[0]*dim
    for i in range(dim):
        new_amp_vec[i]=amp_vec[mapper[i]]
    return new_amp_vec


if __name__ == '__main__':
    print(get_inverse_amp_indices(5))
