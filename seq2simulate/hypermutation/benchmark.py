import hypermutate
import time

def timer(f):
    def decorated(*args, **kwargs):
        start_time = time.time()
        result = f(*args, **kwargs)
        print("func: {}, sec: {}".format(f, time.time() - start_time))
        return result
    return decorated

def get_seq():
    with open('../../data/split/Seq6_Sus', 'r') as f:
        data = f.read()[9:]
    data = data.replace('\n', '')    
    seq = hypermutate.Sequence(data)
    return seq

@timer
def mutate_seq_5(seq):
    seq.mutate_sequence(5)

@timer
def mutate_seq_50(seq):
    seq.mutate_sequence(50)

if __name__ == '__main__':
    seq = get_seq()
    mutate_seq_5(seq)
    seq = get_seq()
    mutate_seq_50(seq)
