from matplotlib import pyplot as plt
from legendre_prng import bruteforce
from numpy import logspace

def main():
    security_bits = 40
    stream_length = 1000000
    keyspace_bits = 22

    print('\ncomputing for variable security bits...')
    xvals = range(30,140,10)
    exec_times = [bruteforce(k, stream_length, keyspace_bits) for k in xvals]
    plot(exec_times, 'log_2(p)', xvals)

    print('\ncomputing for variable stream length...')
    xvals = [int(l) for l in logspace(2.5,7,10)]
    exec_times = [bruteforce(security_bits, k, keyspace_bits) for k in xvals]
    plot(exec_times, 'hint\_bitlength', xvals)

    print('\ncomputing for variable keyspace bits...')
    xvals = range(21, 28)
    exec_times = [bruteforce(security_bits, stream_length, k) for k in xvals]
    plot(exec_times, 'keyspace', xvals)

def plot(data: list, variable: str, x):
    exec_times_v2 = [a for (a,_) in data]
    exec_times_v3 = [b for (_,b) in data]

    p1 = plt.plot(x, exec_times_v2, color='red', marker='o')
    p2 = plt.plot(x, exec_times_v3, color='blue', marker='o')
    plt.legend((p1[0], p2[0]), ('Version 2', 'Version 3'))
    plt.title(r'Execution times depending on \texttt{' + variable + r'}')
    plt.xlabel('Increasing values for $' + variable + r'$')
    plt.ylabel('Time in seconds')
    plt.savefig('figs/inc-'+variable+'-time.svg')
    plt.close()

if __name__ == "__main__":
    main()
