import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import java.util.Random;
import java.util.Properties;
import java.util.Arrays;
import java.util.ArrayList;

public class Individual implements Comparable<Individual>
{
		public double[] values_;
		public double[] sigmas_;
		public double fitness_;
		public int dimensions_;

		public Individual(int dimensions)
		{
			dimensions_ = dimensions;
			values_ = new double[dimensions];
			sigmas_ = new double[dimensions];

		}

		public int getDimensions(){
			return dimensions_;
		}

		public void setValue(int index, double value)
		{
			this.values_[index] = value;
		}
		public double getValue(int index)
		{
			return this.values_[index];
		}

		public void setValues(double[] valueList)
		{
			if (values_.length != valueList.length)
          		throw new IllegalArgumentException ("Cannot combine lists with dissimilar sizes");
          	this.values_ = valueList;
		}

		public double[] getValues()
		{
			return this.values_;
		}

		public void setSigma(int index, double sigma)
		{
			this.sigmas_[index] = sigma;
		}
		public double getSigma(int index)
		{
			return this.sigmas_[index];
		}

		public void setSigmas(double[] sigmaList)
		{
			this.sigmas_ = sigmaList;
		}
		public double[] getSigmas()
		{
			return this.sigmas_;
		}
		

		public void setFitness(double fitness){
			this.fitness_ = fitness;
		}
		public double getFitness()
		{
			return this.fitness_;
		}


		public Individual copyIndividual(){
			Individual newIndividual = new Individual(this.dimensions_);
			newIndividual.setValues(this.values_.clone());
			newIndividual.setSigmas(this.sigmas_.clone());
			return newIndividual;
		}

		@Override
		public int compareTo(final Individual o)
		{
			return Double.compare(this.fitness_, o.fitness_);
		}
}