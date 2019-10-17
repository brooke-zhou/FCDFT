import sys,argparse,ast,re,fnmatch,subprocess,glob,time,scipy
from scipy.optimize import minimize
import numpy as np
import penalty_functions

Niter = 1

def main(argv):
    global Niter

    parser = argparse.ArgumentParser(description=
        'Calculate Fock correction paramters for given elements')

    # arguments
    parser.add_argument('elements', type=str, default='H', nargs='*', 
        help = "Elements in the training set, light -> heavy")
    parser.add_argument('-err', type=str, default='debug', nargs=1, 
        help = "Penalty function type")
    parser.add_argument('-minimizer', type=str, default='BFGS', nargs=1, 
        help = "Optimization algorithm")
    parser.add_argument('-gtol', type=float, default=1e-5, nargs=1, 
        help = "gradient tolerance in BFGS")
    parser.add_argument('-eps', type=float, 
        default=[1e-5,1e-5,1e-5,1e-5,1e-5,1e-5,1e-5], nargs=8, help = "eps in BFGS")
    parser.add_argument('-guess', type=str, default='0', nargs=1, 
        help = "Initial guess option or filename. \
               'zeros': all zero; \
               'step0': FC result of no-Ucor and no-PEC; \
               'kaito': Kaito full; \
               'guess_list.csv': read from 'guess_list.csv' in the same folder")
    
    args = parser.parse_args()

    elements = args.elements
    err = args.err[0]
    minimizer = args.minimizer
    guess = args.guess[0]
    gtol = args.gtol
    eps = args.eps

    if elements[0] == 'H':
        # param = [epsilon_H1s, lambda_H1s]
        nparam = 2
        initial_guess = np.ones(nparam)
        if guess == 'zeros': # initial guess as all zeros
            initial_guess = initial_guess*0.0
        elif guess == 'step0': # initial guess as all zeros
            initial_guess = np.asarray([0.08497394, 0.0])
        elif guess.split('.')[-1] == 'csv': # initial guess is read from file
            initial_guess = np.genfromtxt(guess, delimiter=',')
        else:
            sys.exit('Error: invalid initial guess option!')

        if err == 'square':
            penalty = penalty_functions.H_square_err
        elif err == 'debug':
            penalty = penalty_functions.debug_err
        elif err == 'square_U_only':
            penalty = penalty_functions.H_U_only_square_err
            initial_guess = np.asarray([\
                0.0,0.0,0.0,0.0,0.0,\
                0.0,0.0,0.0,0.0,0.0,\
                0.0,0.0,0.0,0.0,0.0,\
                0.0,0.0,0.0,0.0,0.0\
                ])
        elif err == 'square_U':
            penalty = penalty_functions.H_U_square_err
            initial_guess = np.append(initial_guess,[\
                -0.0323256815,0.0934846152,0.0149072063,-0.1642341001,0.0020415456,\
                -0.0229921778,0.0915471996,-0.0342405310,-0.1634174819,0.0403589439,\
                -0.0072951247,0.0595323756,-0.1226048736,-0.1311303267,0.5012961827,\
                -0.0005398132,0.0107962648,-0.0809719858,0.2699066194,-0.3373832742
                ])
        else:
            sys.exit('Error: invalid penalty function option!')

        if len(elements) >= 2:
            if elements[1] == 'C':
                penalty = penalty_functions.HC_square_err
                # parameters in Ucor not included
                if guess == 'zeros': # initial guess as all zeros
                    if err == 'HC_square_err':
                        initial_guess = np.append(initial_guess, np.ones(5)*0.0)
                    elif err == 'debug':
                        initial_guess = np.append(initial_guess, np.ones(5)*0.0)
                        penalty = penalty_functions.debug_err
                    elif err == 'square':
                        initial_guess = np.append(initial_guess, np.ones(5)*0.0)
                    else:
                        sys.exit('Error: invalid penalty error option!')
                elif guess == 'step0': 
                    if err == 'square':
                        initial_guess = np.asarray([\
                            -0.076846994,-0.911397048,-0.072740548,-0.137419801,\
                            0.006646586,0.242750117,-0.246016173])
                    elif err == 'debug':
                        initial_guess = np.asarray([\
                            -0.08497394,-0.053817519,-0.005272491,-0.086407013,\
                            0.0, 0.0, 0.0, 0.0])
                        penalty = penalty_functions.debug_err
                    else:
                        sys.exit('Error: invalid penalty function option!')
                elif guess.split('.')[-1] == 'csv': # initial guess is read from file
                    initial_guess = np.genfromtxt(guess, delimiter=',')
                else:
                    sys.exit('Error: invalid initial guess option!')
            else:
                sys.exit('Error: invalid element!')

        if len(elements) >= 3:
            if elements[2] == 'O':
                initial_guess = np.append(initial_guess, np.ones(0)*guess)
                if err == 'debug':
                    penalty = penalty_functions.debug_err
                elif err == 'square':
                    penalty = penalty_functions.HCO_square_err
                else:
                    sys.exit('Error: invalid penalty error option!')
            else:
                sys.exit('Error: invalid element!')
    else:
        sys.exit('Error: invalid element!')

    # callback function to print convergence information
    def callbackH(Xi):
        global Niter
        fval = penalty(Xi)
        print( ' - penalty function at this step F(x)    = ' + str(fval) )
        f.write( '{0:4d}     {1: 3.10f}     {2: 3.10f}     {3: 3.10f}\n'.\
            format(Niter, Xi[0], Xi[1], fval) )
        Niter += 1
    def callbackHC(Xi):
        global Niter
        fval = penalty(Xi)
        print( ' - penalty function at this step F(x) = ' + str(fval) )
        # f.write( '{0:4d}     {1: 3.10f}\n'.format(Niter, fval) ) # save only penalty function to file
        f.write( '{0:4d}     {1: 3.10f}     {2: 3.10f}     {3: 3.10f}     {4: 3.10f}     {5: 3.10f}     {6: 3.10f}     {7: 3.10f}     {8: 3.10f}     \n'.\
                            format(Niter, Xi[0], Xi[1], Xi[2], Xi[3], Xi[4], Xi[5], Xi[6], fval) )
        Niter += 1

    # BFGS
    # jac: Function returns gradient or not.
    # tol : Optional Tolerance for termination. If not None, will be the default value for gtol
    # maxiter : Maximum number of iteration. if maxiter is None: maxiter = len(x0) * 200
    # disp : Set to True to print convergence messages
        # [Current function value,Iterations,Function evaluations,Gradient evaluations]
    # gtol : Gradient norm must be less than gtol before successful termination. 
    # norm : Order of norm.
    # eps : Step size for numerically calculating the gradient. float or ndarray
    # return_all : Return a list of x at each iteration if True.
    # For more information, see scipy docs

    if guess.split('.')[-1] == 'csv': # list of initial guesses 
        n_guess = 0
        for iguess in initial_guess:
            n_guess += 1
            filename = 'fc_init-no' + str(n_guess) + '.txt'
            # with open(filename,'w') as f:
            #     f.write( 'Training Set: ' + ", ".join(str(x) for x in elements) + '\n' )
            #     f.write( 'Error function type: ' + err + '\n' )
            #     f.write( 'Minimizer type: ' + minimizer + '\n' )
            #     f.write( 'Initial guess option is: ' + str(guess) + '\n' )
            #     Niter = 1
            #     start_time = time.time()
            #     if len(elements) == 1:
            #         f.write( '-------------------------------------------------------------------------------------\n' )
            #         f.write(    '{0:4s}     {1:13s}     {2:13s}     {3:13s}     {4:13s}     {5:13s}\n'.\
            #                             format('Iter', ' epsilon', ' alpha', ' beta', ' gamma', ' error') )
            #         f.write( '{0:4d}     {1: 3.10f}     {2: 3.10f}     {3: 3.10f}     {4: 3.10f}     {5: 3.10f}\n'.\
            #                             format(0, iguess[0], iguess[1], iguess[2], 2.0, penalty(iguess)) )
            #         q = scipy.optimize.minimize(penalty, iguess, args=(), method=minimizer, jac=None, tol=None, callback=callbackH, \
            #             options={'disp': False, 'gtol': gtol, 'eps': eps, 'return_all': True, 'maxiter': 100, 'norm': float('inf')})
            #         f.write( '-------------------------------------------------------------------------------------\n' )
            #     elif len(elements) is 2 and err is not 'HC_diag_square_err':
            #         # f.write( '--------------------\n' )
            #         f.write( '----------------------------------------------------------------------------------------------------------------------------------------------------\n' )
            #         f.write(    '{0:4s}     {1:13s}\n'.format('Iter', ' error') )
            #         f.write( '{0:4d}     {1: 3.10f}\n'.format(0, penalty(initial_guess)) )
            #         q = scipy.optimize.minimize(penalty, iguess, args=(), method=minimizer, jac=None, tol=None, callback=callbackHC, \
            #             options={'disp': False, 'gtol': gtol, 'eps': eps, 'return_all': False, 'maxiter': 100, 'norm': float('inf')})
            #         f.write( '----------------------------------------------------------------------------------------------------------------------------------------------------\n' )
            #         # f.write( '--------------------\n' )
            #     elif len(elements) == 2 and err is 'HC_diag_square_err':
            #         f.write( '------------------------------------------------------------------------------------------------------------\n' )
            #         f.write(    '{0:4s}     {1:13s}     {2:13s}     {3:13s}     {4:13s}     {5:13s}\n'.\
            #                            format('Iter', ' e_H1s', ' e_C1s', ' e_C2s', ' e_C2p', ' error') )
            #         f.write( '{0:4d}     {1: 3.10f}     {2: 3.10f}     {3: 3.10f}     {4: 3.10f}     {5: 3.10f}\n'.\
            #                            format(0, iguess[0], iguess[1], iguess[2], iguess[3], penalty(initial_guess)) )
            #         q = scipy.optimize.minimize(penalty, iguess, args=(), method=minimizer, jac=None, tol=None, callback=callbackHC_diag, \
            #             options={'disp': False, 'gtol': gtol, 'eps': eps, 'return_all': False, 'maxiter': 100, 'norm': float('inf')})
            #         f.write( '------------------------------------------------------------------------------------------------------------\n' )
                    
            #     elapsed_time = time.time() - start_time
            #     f.write( 'Total optimization time : ' + str(elapsed_time) + '\n' )
            #     f.write( 'Convergence status is : ' + str(q.success) + '\n' )
            #     f.write( q.message + '\n' )
            #     f.write( 'Total number of function evaluations: ' + str(q.nfev) + '\n' )
            #     f.write( 'Total number of derivative evaluations: ' + str(q.njev) + '\n' )
            #     f.write( 'The final parameters are: ' + str(q.x) + '\n' )
            #     np.savetxt('allgrads' + '_init-no' + str(n_guess) +'.csv', q.allgrads, delimiter=',' )
            with open(filename,'w') as f:
                f.write( 'Training Set: ' + '[' + ", ".join(str(x) for x in elements) + ']' + '\n' )
                f.write( 'Error function type: ' + err + '\n' )
                f.write( 'Minimizer type: ' + minimizer + '\n' )
                f.write( 'Initial guess option is: ' + '[' + ", ".join(str(x) for x in iguess) + ']' + '\n' )
                f.write( 'Gradient threshold is gtol = ' + str(gtol) + '\n' )
                f.write( 'Gradient stepsize is eps = ' + ", ".join(str(x) for x in eps) + '\n' )
                Niter = 1
                start_time = time.time()
                if len(elements) == 1:
                    f.write( '-------------------------------------------------------------------------------------\n' )
                    f.write(    '{0:4s}     {1:13s}     {2:13s}     {3:13s}\n'.\
                                        format('Iter', ' epsilon_H1s', ' lambda_H1s', ' error') )
                    f.write( '{0:4d}     {1: 3.10f}     {2: 3.10f}     {3: 3.10f}\n'.\
                                        format(0, initial_guess[0], initial_guess[1], initial_guess[2], 2.0, penalty(initial_guess)) )
                    q = scipy.optimize.minimize(penalty, initial_guess, args=(), method=minimizer, jac=None, tol=None, callback=callbackH, \
                        options={'disp': False, 'gtol': gtol, 'eps': eps, 'return_all': True, 'maxiter': 100, 'norm': float('inf')})
                    f.write( '-------------------------------------------------------------------------------------\n' )
                elif len(elements) == 2:
                    f.write( '-----------------------------------------------------------------------------------------------------------------------------------\n' )
                        # f.write( '--------------------\n' ) 
    
                    # f.write(    '{0:4s}     {1:13s}\n'.format('Iter', ' error') ) # print only penalty function
                    # f.write( '{0:4d}     {1: 3.10f}\n'.format(0, penalty(initial_guess)) ) # print only penalty function
                    f.write(    '{0:4s}     {1:13s}     {2:13s}     {3:13s}     {4:13s}     {5:13s}     {6:13s}     {7:13s}     {8:13s}\n'.\
                                    format('Iter',' epsilon_H1s', ' epsilon_C1s', ' epsilon_C2s', ' epsilon_C2p', \
                                                    ' lambda_H1s', ' lambda_C2s', ' lambda_C2p', ' error') ) 
                    f.write( '{0:4d}     {1: 3.10f}     {2: 3.10f}     {3: 3.10f}     {4: 3.10f}     {5: 3.10f}     {6: 3.10f}     {7: 3.10f}     {8: 3.10f}\n'.\
                                        format(0, iguess[0], iguess[1], iguess[2], iguess[3], \
                                                            iguess[4], iguess[5], iguess[6], penalty(iguess)) ) # print params and penalty function
    
                    q = scipy.optimize.minimize(penalty, iguess, args=(), method=minimizer, jac=None, tol=None, callback=callbackHC, \
                        options={'disp': False, 'gtol': gtol, 'eps': eps, 'return_all': True, 'maxiter': 100, 'norm': float('inf')})
                    f.write( '-----------------------------------------------------------------------------------------------------------------------------------\n' )
                    # f.write( '--------------------\n' )
    
                elapsed_time = time.time() - start_time
                f.write( 'Total optimization time : ' + str(elapsed_time) + '\n' )
                f.write( 'Convergence status is : ' + str(q.success) + '\n' )
                f.write( q.message + '\n' )
                f.write( 'Total number of function evaluations: ' + str(q.nfev) + '\n' )
                f.write( 'Total number of derivative evaluations: ' + str(q.njev) + '\n' )
                f.write( 'The final parameters are: ' + str(q.x) + '\n' )
                np.savetxt('g' + '_init-no' + str(n_guess) +'.csv', q.allgrads, delimiter=',' )

    else: # single-value initial_guess
        filename = 'fc_init-' + guess + '.txt'
        with open(filename,'w') as f:
            f.write( 'Training Set: ' + ", ".join(str(x) for x in elements) + '\n' )
            f.write( 'Error function type: ' + err + '\n' )
            f.write( 'Minimizer type: ' + minimizer + '\n' )
            f.write( 'Initial guess option is: ' + str(guess) + '\n' )
            f.write( 'Gradient threshold is gtol = ' + str(gtol) + '\n' )
            f.write( 'Gradient stepsize is eps = ' + ", ".join(str(x) for x in eps) + '\n' )
            Niter = 1
            start_time = time.time()
            if len(elements) == 1:
                f.write( '-------------------------------------------------------------------------------------\n' )
                f.write(    '{0:4s}     {1:13s}     {2:13s}     {3:13s}\n'.\
                                    format('Iter', ' epsilon_H1s', ' lambda_H1s', ' error') )
                f.write( '{0:4d}     {1: 3.10f}     {2: 3.10f}     {3: 3.10f}\n'.\
                                    format(0, initial_guess[0], initial_guess[1], initial_guess[2], 2.0, penalty(initial_guess)) )
                q = scipy.optimize.minimize(penalty, initial_guess, args=(), method=minimizer, jac=None, tol=None, callback=callbackH, \
                    options={'disp': False, 'gtol': gtol, 'eps': eps, 'return_all': True, 'maxiter': 100, 'norm': float('inf')})
                f.write( '-------------------------------------------------------------------------------------\n' )
            elif len(elements) == 2:
                f.write( '-----------------------------------------------------------------------------------------------------------------------------------\n' )
                # f.write( '--------------------\n' ) 

                # f.write(    '{0:4s}     {1:13s}\n'.format('Iter', ' error') ) # print only penalty function
                # f.write( '{0:4d}     {1: 3.10f}\n'.format(0, penalty(initial_guess)) ) # print only penalty function
                f.write(    '{0:4s}     {1:13s}     {2:13s}     {3:13s}     {4:13s}     {5:13s}     {6:13s}     {7:13s}     {8:13s}\n'.\
                                format('Iter',' epsilon_H1s', ' epsilon_C1s', ' epsilon_C2s', ' epsilon_C2p', \
                                                ' lambda_H1s', ' lambda_C2s', ' lambda_C2p', ' error') ) 
                f.write( '{0:4d}     {1: 3.10f}     {2: 3.10f}     {3: 3.10f}     {4: 3.10f}     {5: 3.10f}     {6: 3.10f}     {7: 3.10f}     {8: 3.10f}\n'.\
                                    format(0, initial_guess[0], initial_guess[1], initial_guess[2], initial_guess[3], \
                                                        initial_guess[4], initial_guess[5], initial_guess[6], penalty(initial_guess)) ) # print params and penalty function

                q = scipy.optimize.minimize(penalty, initial_guess, args=(), method=minimizer, jac=None, tol=None, callback=callbackHC, \
                    options={'disp': False, 'gtol': gtol, 'eps': eps, 'return_all': True, 'maxiter': 100, 'norm': float('inf')})
                f.write( '-----------------------------------------------------------------------------------------------------------------------------------\n' )
                # f.write( '--------------------\n' )

            elapsed_time = time.time() - start_time
            f.write( 'Total optimization time : ' + str(elapsed_time) + '\n' )
            f.write( 'Convergence status is : ' + str(q.success) + '\n' )
            f.write( q.message + '\n' )
            f.write( 'Total number of function evaluations: ' + str(q.nfev) + '\n' )
            f.write( 'Total number of derivative evaluations: ' + str(q.njev) + '\n' )
            f.write( 'The final parameters are: ' + str(q.x) + '\n' )

            np.savetxt( 'g.csv', q.allgrads, delimiter=',' )

# Allows the script to run as an exe
if __name__ == "__main__":
     main(sys.argv[1:])


