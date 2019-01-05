def play():
    while True:
        name = raw_input('Please type your name here \n')
        if name == 'Yiran':
            print('Yeahhhhh!!!!! Try typing something else.')
            continue
        elif name == 'break':
            break
        else:
            print('That wasn\'t your name!')
