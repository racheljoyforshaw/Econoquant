
# coding: utf-8

# Simulates the private learning element of Cripps & Thomas (2014) "Strategic Experimentation in Queues".  
# Paper available here: http://www.ed.ac.uk/polopoly_fs/1.150007!/fileManager/queues2402131.pdf


#IMPORT SOME USEFUL PACKAGES
import numpy as np
import math
import matplotlib.pyplot as plt



#SET UP A CLASS OF QUEUE MEMBERS
class Queuer:
    """ Class of queue members """
    
    
    def __init__(self, queue_length, tau):
        # assigns (improper) prior about state of server is drawn from the uniform distribution 
        self.mu_prior = 0.99 #np.random.uniform() #0.99
        # arrival time
        self.tau = tau
        # queue position
        self.queue_position = queue_length + 1.0
        # intitial balk point
        self.mu_cutoff = (1.0-delta)/delta*alpha*(psi**(self.queue_position-1)*w - 1.0)
        
        
    def update_prior(self, tau, service):
        if service == 1:
            #observe service
            self.mu_posterior = 1.0
        else:
            #has observed service while queuing
            try:
                if self.mu_posterior == 1.0:
                    pass
                else:
                #has not yet observed service but has has an update before
                    self.mu_posterior = (self.mu_prior*(1.0-alpha)**(tau-self.tau))/                                          (self.mu_prior*(1.0-alpha)**(tau-self.tau) + 1.0 - self.mu_prior)
            except:
                # has not observed service, has not had an update before
                self.mu_posterior = (self.mu_prior*(1.0-alpha)**(tau-self.tau))/                                      (self.mu_prior*(1.0-alpha)**(tau-self.tau) + 1.0 - self.mu_prior)
    
    def update_cutoff(self,removals):
        self.queue_position = self.queue_position - removals
        self.mu_cutoff = (1-delta)/(delta*alpha*psi**(self.queue_position-1.0)*w - 1.0)



def queueBuilder(delta, alpha, balk_payoff, w, k, T, quiet):
    """ Builds a queue of length T given inputs:
    delta = discount rate
    alpha = service rate
    balk_payoff = outside option
    w = payoff from service
    k = possible number of people served per service period
    quiet = 1 surpresses the log  
    """
    # preliminaries...
    psi = alpha/(1.0-delta*(1.0-alpha)) # congestion cost
    Queue = []
    queue_length = 0
    max_queue_length = 0
    queue_sim = []
    total_balkers = 0
    
    # lets build a queue...
    for tau in xrange(1,T):
        if quiet ==0:
            print 'tau = ' + str(tau) + ':  ' + str(queue_length) + ' in queue' 
        
        # SERVICE
        # z is number of trials to get one success drawn once from geometric distribution
        z =  np.random.geometric(alpha,1)
        if z == 1:
            # if number of trials for sucess = 1, enter k into queue
            y = 1.0 * k # Is this right for k>1 ???!!!
        else:
            y = 0.0 
   
        # take y people out of the queue - they have been served
        if y >0:
            if quiet == 0:
                print '         -' + str(int(y)) + ' served'
        else:
            if quiet == 0:
                print '          ' + str(int(y)) + ' served'
        if queue_length == 0:
            #no-one in queue - do nothing  
            pass
    
        # no service
        elif y == 0.0:
            for i in range(queue_length):
                #update priors and cutoffs
                Queue[i].update_prior(tau, service = 0)
                Queue[i].update_cutoff(0.0)
            
        #more or equal removals than people in queue
        elif queue_length <= int(y):
            #get rid of all queuers
            for i in range(queue_length):
                Queue.remove(Queue[i])
            
            queue_length = 0

            
        #fewer removals than people in queue
        else:
            queue_length = queue_length - int(y)
            # remove those who have been served
            for i in range(int(y)):
                Queue.remove(Queue[i])
            #update the priors and cutoffs of those remaining in the queue    
            for i in range(queue_length):
                Queue[i].update_prior(tau,service = 1)
                Queue[i].update_cutoff(y)

        # EXIT
        balkers = 0

    
        # DECIDE WHETHER OR NOT TO BALK
        for i in range(queue_length):
            if quiet == 0:
                print str(int(Queue[i].queue_position)) + ' in queue:'
            
            try:
                if quiet ==0:
                    print 'prior = ' + str(Queue[i].mu_prior)
                    print 'posterior = ' + str(Queue[i].mu_posterior)
                    print 'cutoff = ' + str(Queue[i].mu_cutoff)
            
                if Queue[i].mu_posterior <= Queue[i].mu_cutoff:
                    #if posterior falls below cutoff - balk
                    Queue.remove(Queue[i])
                    balkers = balkers + 1
                    queue_length = queue_length - 1
                    #update the cutoffs of those behind the balker in the queue    
                    for k in range(i+1,queue_length):
                        Queue[k].update_cutoff(1)


            except:
                # they don't have a posterior yet
                pass
        
    
        total_balkers = total_balkers + balkers
        
        if balkers > 0:
            if quiet ==0:
                print '         -' + str(balkers) + ' balked'
            
        else:
            if quiet == 0:
                print '          ' + str(balkers) + ' balked'

        # ARRIVAL
        x = Queuer(queue_length, tau)
        # balk before joining if outside payoff higher
        if psi**(x.queue_position+1)*delta*w > balk_payoff:
            #join
            Queue.append(x)
            queue_length = queue_length + 1
            if quiet ==0:
                print '         +1 arrival'
        else:
            #balk
            if quiet==0:
                print '         +0 arrival'
    
        #keep track of the maximum queue length
        if queue_length> max_queue_length:
            max_queue_length = queue_length
    
        #make a nice array of queue length at given tau
        queue_sim.append(queue_length)

    #spit these out please
    return max_queue_length,queue_sim
        

# PARAMETERS
delta = 0.8 # discount rate
alpha =  0.5 # service rate
balk_payoff = 1 # payoff from balking
w = 10.0 # payoff from service
k = 1.0 # number served
psi = alpha/(1.0-delta*(1.0-alpha)) # congestion cost
T = 100 # how many periods of queueing?



[max_queue_length, queue_sim] = queueBuilder(delta, alpha, balk_payoff, w, k, T, quiet=1)
    


# Plot the queue - where the queue levels out is M
plt.plot(queue_sim)
plt.ylim(0,max_queue_length+1)
plt.xlim(0, T)
plt.show



# Let's build lots of queues! (because we have specified a stochastic mu)
iterations = 20
Q = np.empty((iterations,T-1,2)) # how many queues?
for j in range(iterations):
    # lets build a queue...
    [max_queue_length, queue_sim] = queueBuilder(delta, alpha, balk_payoff, w, k, T, quiet=1)
    Q[j,:,1] = queue_sim



#...and plot them
for i in range(iterations):
    plt.plot(Q[i,:,1])
    plt.ylim(0,max_queue_length+1)
    plt.xlim(0, T)
    plt.show

