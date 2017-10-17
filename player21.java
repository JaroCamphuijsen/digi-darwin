import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import java.util.Random;
import java.util.Properties;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.ArrayList;
import java.lang.NullPointerException;




public class player21 implements ContestSubmission
{
	Random rnd_;
	ContestEvaluation evaluation_;
    private int evaluations_limit_;
    private int nParents_;
    private int popSize_;
    private int nChildren_;
    private double pSelFactor_;
    private boolean recombSwitch_;
    private int nMutations_;
    private double sigBound_;
	
	public player21()
	{
		rnd_ = new Random();
	}
	
	public void setSeed(long seed)
	{
		// Set seed of algortihms random process
		rnd_.setSeed(seed);
	}

	public void setEvaluation(ContestEvaluation evaluation)
	{
		// Set evaluation problem used in the run
		evaluation_ = evaluation;
		
		// Get evaluation properties
		Properties props = evaluation.getProperties();
        // Get evaluation limit
        evaluations_limit_ = Integer.parseInt(props.getProperty("Evaluations"));
		// Property keys depend on specific evaluation
		// E.g. double param = Double.parseDouble(props.getProperty("property_name"));
        boolean isMultimodal = Boolean.parseBoolean(props.getProperty("Multimodal"));
        boolean hasStructure = Boolean.parseBoolean(props.getProperty("Regular"));
        boolean isSeparable = Boolean.parseBoolean(props.getProperty("Separable"));

		// Do sth with property values, e.g. specify relevant settings of your algorithm

        if(!isMultimodal){
        	// Changing this will affect Sphere
        	// Changing this will affect function 1

        	this.popSize_ = 20;
            this.nParents_ = this.popSize_*2;
            this.nChildren_ = this.nParents_;
        	this.pSelFactor_ = .3;
        	this.pSelFactor_ = 1.0;

        	this.recombSwitch_= true;
        	this.nMutations_ = 1;
        	this.sigBound_ = 0.1;
        }
        else if(hasStructure){
        	// Changing this will affect Rastrigin
        	// Changing this will affect function 3

        	this.popSize_ = 250;
            this.nParents_ = this.popSize_*2;
            this.nChildren_ = (int)Math.ceil(this.nParents_*1);
        	this.pSelFactor_ = 1.0;
        	this.recombSwitch_= true;
        	this.nMutations_ = 2;
        	this.sigBound_ = 0.5;

        }
        else {
        	// Changing this will affect Fletcher Powell
        	// Changing this will affect function 2
        	this.popSize_ = 200; 
            this.nParents_ = this.popSize_ * 2 ;	
            this.nChildren_ = this.nParents_;
            this.pSelFactor_ = .02; 
        	this.recombSwitch_= true;
        	this.nMutations_ = 1;
        	this.sigBound_ = .1; 
        }
        if(isSeparable){
        	this.nMutations_ = 1;
        }
    }
    
	public void run()
	{
		// Run your algorithm here
        
        int evals = 0;
        int dimensions = 10;
        int popSize = this.popSize_;
        int nParents = this.nParents_;
        int nChildren = this.nChildren_;
        int nMutations = this.nMutations_;
        double a = -5.0;
        double b = 5.0;
        double pSelFactor = this.pSelFactor_; // make smaller (than one) to increase variation in chosen parents
		double sigBound = this.sigBound_; // lower bound on sigma in new mutants
		double pSelFactorInit = pSelFactor;
		double sigBoundInit = sigBound;

       	// init population
		List<Individual> population = new ArrayList<Individual>();
		List<Individual> matingPool = new ArrayList<Individual>();
		List<Individual> offspring = new ArrayList<Individual>();
        for (int i=0; i<popSize; i++){
        	population.add(new Individual(dimensions));
        	for (int j=0; j<dimensions; j++){
        		double val = rnd_.nextDouble() * (b - a) + a;
	        	population.get(i).setValue(j, val);
	        	double sig = rnd_.nextGaussian() * 0.01;
	        	population.get(i).setSigma(j, sig);
        	}
        	double[] values = population.get(i).getValues();

        	double fitness = (double) evaluation_.evaluate(values);
        	evals++;
        	population.get(i).setFitness(fitness);
        }

        //Computing the Cummulative Probability Distribution, if we take a ranking based probability we can predefine the whole distribution
        double[] cpb = new double[popSize];
			double probSum = 0;
			for (int i=0; i<popSize;i++){
				probSum = probSum + (1 - Math.exp(i * pSelFactor));
				cpb[i] = probSum;
			}

			//Normalize probability distribution
			for (int i=0; i<popSize;i++){
				cpb[i] = cpb[i]/probSum;

			}



		    // Sort population on their fitness
        	Collections.sort(population);



        // Start evaluation cycle
        double wheelStep;
        double evalCount = 0;
        outerloop:
        while(evals<evaluations_limit_){
    	// Set parameters that are dependent on time
        	double evalFrac = (double) evals/evaluations_limit_;
        	sigBound = sigBoundInit*Math.pow((1-evalFrac),4);
        	pSelFactor = pSelFactorInit*Math.pow((evalFrac),3);
        	// System.out.println(pSelFactor);

        	if (evalFrac>evalCount){
        		System.out.print("|");
	        	evalCount = evalCount + 0.05; 
			}


			//Computing the Cummulative Probability Distribution with time dependent pSelFactor
        	cpb = new double[popSize];
			probSum = 0;
			for (int i=0; i<popSize;i++){
				probSum = probSum + (1 - Math.exp(i * pSelFactor));
				cpb[i] = probSum;
			}

			//Normalize probability distribution
			for (int i=0; i<popSize;i++){
				cpb[i] = cpb[i]/probSum;
			}



		    // Sort population on their fitness
        	Collections.sort(population);

			// Stochastic Universal Sampling
			int cpbIndex = 0;
			int currentMember = 0;
			matingPool = new ArrayList<Individual>();
			double[] parentCounter = new double[popSize];
			for (int n=0; n<popSize;n++){
				parentCounter[n] = 0;
			}
			for (int n=0; n<nParents;n++){
				matingPool.add(new Individual(dimensions));
			}

			wheelStep = rnd_.nextDouble() * (double) 1/nParents;


			while (currentMember<nParents){
				innerloop:
			 	while (wheelStep < cpb[cpbIndex]){

			 		matingPool.set(currentMember, population.get(cpbIndex).copyIndividual());
			 		wheelStep += (double) 1/nParents;
			 		parentCounter[cpbIndex] += 1;
			 		currentMember += 1;
			 		if (currentMember >= nParents){
			 			break innerloop;
			 		}
			 	}
			 	cpbIndex += 1;
			}


            // Reproduce --> make offspring

            offspring = new ArrayList<Individual>();
            for (int i = 0; i < nChildren; i++){
            	Individual child = new Individual(dimensions);

				//Apply recombination operator to obtain new children
				if (this.recombSwitch_){
            		child = recombination(matingPool.get(rnd_.nextInt(matingPool.size())), matingPool.get(rnd_.nextInt(matingPool.size())));
            	} else{
            		child = matingPool.get(rnd_.nextInt(matingPool.size()));
            	}
            	// Apply mutation operator on the new children
				for (int m = 0; m<nMutations;m++){
					mutation(child, sigBound);
				}
				// Evaluate the fitness of the child and add to offspring
				try{
					child.setFitness((double) evaluation_.evaluate(child.getValues()));
					evals++;
				}
				catch(Exception e){
            		System.out.println("MAX EVALUATIONS REACHED!");
            		break outerloop;
				}
				

				offspring.add(child);
        	}
        	

			// Select survivors
            Collections.sort(offspring);
            
     
            population = new ArrayList<Individual>();
            for (int i = 0; i < popSize; i++){
            	population.add(offspring.get(offspring.size() - (i + 1)));
            }
            Collections.reverse(population);

            // calculate popmean
            double popSum = 0;
             for (int i = 0; i < popSize; i++){
             	popSum += population.get(i).getFitness();
            }

        }
	}

	//mutation operator self adapting multiple sigma
	public void  mutation(Individual mutant, double sigBound)
	{	
		double tau = 1/ Math.sqrt(2*mutant.getDimensions());
		double tau_1 = 1/ Math.sqrt(2*Math.sqrt(mutant.getDimensions()));

		int index = rnd_.nextInt(mutant.getDimensions());
		double sigma = mutant.getSigma(index) * Math.exp(tau_1 * rnd_.nextGaussian() + tau* rnd_.nextGaussian());
		sigma = Math.max(sigma, sigBound);
		mutant.setSigma(index, sigma);
		double value = mutant.getValue(index) + sigma * rnd_.nextGaussian();
		mutant.setValue(index, value);

	}

	//uniform arithmetic recombination
	public Individual recombination(Individual parent1, Individual parent2)
	{
		Individual child = parent1.copyIndividual();
		for (int i = 0; i<child.getDimensions(); i++){
			double value = .5 * parent1.getValue(i) + .5 * parent2.getValue(i);
			child.setValue(i, value);
			double sigma = .5 * parent1.getSigma(i) + .5 * parent2.getSigma(i);
			child.setSigma(i, sigma);
		}

		return child;
	}

	
	// Blend Crossover
		public Individual recombinationBLX(Individual parent1, Individual parent2)
	{
		Individual child = parent1.copyIndividual();
		for (int i = 0; i<child.getDimensions(); i++){
			double alpha = .5;
			double u = rnd_.nextDouble();
			double gamma = (1 - 2*alpha) * u - alpha;
			double value = 0;
			double sigma = 0;
			if (parent1.getValue(i)<parent2.getValue(i)){
				value = (1 - gamma) * parent1.getValue(i) + gamma * parent2.getValue(i);
			}
			else{
				value = (1 - gamma) * parent2.getValue(i) + gamma * parent1.getValue(i);

			}
			child.setValue(i, value);

			if (parent1.getSigma(i)<parent2.getSigma(i)){
				sigma = (1 - gamma) * parent1.getSigma(i) + gamma * parent2.getSigma(i);
			}
			else{
				sigma = (1 - gamma) * parent2.getSigma(i) + gamma * parent1.getSigma(i);

			}
			child.setSigma(i, sigma);

		}
		
		return child;
	}
	
}
