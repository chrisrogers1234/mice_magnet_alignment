import json

import fit_tm 

def generate_test_data():
    fout = open('test_data.json', 'w')
    fitter = fit_tm.FitTM()
    test_alignment = {
          "dx":4.0,
          "dy":3.0,
          "dz":2.0,
          "dxp":0.001,
          "dyp":0.002,
          "scale_factor":93.0
        }

    test_data = []
    print "Generating test data..."
    for tof12 in [32., 31.5, 31.0]:
        my_dict = {
            'transfer_matrix':fitter.get_tm(test_alignment, tof12),
            'error_matrix':[[1. for i in range(5)] for j in range(4)], 
            'tof12':tof12,
        }
        test_data.append(my_dict)
        print "TOF12", tof12
        print "Equivalent pz", fitter.tracking.pz
        print "TM:"
        for row in my_dict['transfer_matrix']:
            print "   ", row
    print >> fout, json.dumps(test_data, indent=2)

def test_fitter():
    """
    Initialise the fitter and run it; check that it converges on reasonable
    values.
    """
    fin = open('test_data.json')
    test_data = json.loads(fin.read())

    fitter = fit_tm.FitTM()
    fitter.fit(test_data)
        
def main():
    generate_test_data()
    test_fitter()

if __name__ == "__main__":
    main()
    print "Finished"

