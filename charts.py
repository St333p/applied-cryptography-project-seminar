from matplotlib import pyplot as plt
from legendre_prng import bruteforce
from numpy import logspace

TIME = 0
SYMBOLS = 1

def main():
    security_bits = 40
    stream_length = 1000000
    keyspace_bits = 22

    print('\ncomputing for variable security bits...')
    xvals = range(30,140,10)
    data = [bruteforce(k, stream_length, keyspace_bits) for k in xvals]
    data_time = [d[TIME] for d in data]
    plot(data_time, 'log_2(p)', xvals, TIME)
    data_symbols = [d[SYMBOLS] for d in data]
    plot(data_symbols, 'log_2(p)', xvals, SYMBOLS)

    print('\ncomputing for variable stream length...')
    xvals = [int(l) for l in logspace(2.2,7,10)]
    data = [bruteforce(security_bits, k, keyspace_bits) for k in xvals]
    data_time = [d[TIME] for d in data]
    plot(data_time, 'hintbitlength', xvals, TIME, True)
    data_symbols = [d[SYMBOLS] for d in data]
    plot(data_symbols, 'hintbitlength', xvals, SYMBOLS, True)

    print('\ncomputing for variable keyspace bits...')
    xvals = range(18, 29)
    data = [bruteforce(security_bits, stream_length, k) for k in xvals]
    data_time = [d[TIME] for d in data]
    plot(data_time, 'keyspace', xvals, TIME)
    data_symbols = [d[SYMBOLS] for d in data]
    plot(data_symbols, 'keyspace', xvals, SYMBOLS)

def plot(data: list, variable: str, x, type_y, xlog=False):
    data_v2 = [a for (a,_) in data]
    data_v3 = [b for (_,b) in data]
    
    if type_y == TIME:
        y_str = "execution_time"
    elif type_y == SYMBOLS:
        y_str = "legendre_symbols"
    else:
        raise ValueError("unknown value for type_y")
    
    if xlog:
        plt.xscale('log')
    p1 = plt.plot(x, data_v2, color='red', marker='o')
    p2 = plt.plot(x, data_v3, color='blue', marker='o')
    plt.legend((p1[0], p2[0]), ('Version 2', 'Version 3'))
    plt.title(y_str + r' depending on $' + variable + r'$')
    plt.xlabel('Increasing values for $' + variable + r'$')
    plt.ylabel(y_str)
    plt.savefig('figs/inc-' + variable + '-' + y_str + '.svg')
    plt.close()

if __name__ == "__main__":
    main()
