// inspired from SnAPTreeLikelihood.java of Bryant and Bouckaerts
// modified by CE RAbier for SnappNet 
//computes the network likelihood
// at this time,  ascSiteCountInput and useLogLikelihoodCorrection have not been handled... 


package snappNetForSimSnappNet.core;


import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import javax.xml.bind.JAXBException;

import org.xml.sax.SAXException;

import com.google.common.collect.Multiset; 

import beast.app.BeastMCMC;
import beast.core.Citation;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.TreeInterface;
//import snap.Data;
//import snap.NodeData;
//import snap.likelihood.SnAPLikelihoodCore;
import beast.util.Randomizer;

 
public class SnappNetNetworkLikelihood extends Distribution {
	
	public Input<SnapData> m_pDataInput = new Input<SnapData>("data", "set of alignments");
//	public Input<Tree> m_pTree = new Input<Tree>("tree", "tree with phylogenetic relations");
	public final Input<Network> speciesNetworkInput =
            new Input<>("speciesNetwork", "The species network.", Validate.REQUIRED);
	
	public Input<RealParameter> m_pU = new Input<RealParameter>("mutationRateU", "Instantaneous rate of mutating from the 0 allele to the 1 alelle");
	public Input<RealParameter> m_pV = new Input<RealParameter>("mutationRateV", "Instantaneous rate of mutating from the 1 allele to the 0 alelle");
	public Input<RealParameter> m_pCoalescenceRate = new Input<RealParameter>("coalescenceRate", "population size parameter with one value for each node in the tree");
	public Input<Boolean> m_usenNonPolymorphic = new Input<Boolean>("non-polymorphic", 
			"Check box only if constant sites have been left in the data and are to be included in the likelihood calculation. " +
			"Leave unchecked if all but the variable sites have been removed.",
			//"Whether to use non-polymorphic data in the sequences. " +
			//"If true (default), constant-sites in the data will be used as part of the likelihood calculation. " +
			//"If false , constant sites will be removed from the sequence data and a normalization factor is " +
			//"calculated for the likelihood.", 
			true); 
	
	public Input<IntegerParameter> ascSiteCountInput = new Input<IntegerParameter>("ascSiteCount", "Counts for number of ascertained sites");
	public Input<Boolean> useLogLikelihoodCorrection = new Input<Boolean>("useLogLikelihoodCorrection", "use correction of log likelihood for the purpose of calculating " +
			"Bayes factors for different species assignments. There is (almost) no computational cost involved for the MCMC chain, but the log likelihood " +
			"might be reported as positive number with this correction since the likelihood is not a proper likelihood any more.", false);
	
	
	// At this time, i have not handled ascSiteCount and not useLogLikelihoodCorrection 
	
	
	/** some variable for shadowing inputs **/
	boolean m_bUsenNonPolymorphic;  //true if we use nonPolymorphic sites
	//boolean m_bMutationOnlyAtRoot;
	//boolean m_bHasDominantMarkers;
	
	int numPatterns=0; //will contain the number of patterns (all the sites are summed up in patterns)
	double [] patternProbs; //will contain the likelihood for all sites summed up in patterns
	
	SnapData m_pData; //the data
	Network speciesNetwork=null;
	
	double m_fP0 = 0.0, m_fP1 = 0.0;
	double ascLogP = Double.NaN;
	
	// Correction so that the returned value is a likelihood instead
	// of a sufficient statistic for the likelihood
	double m_fLogLikelihoodCorrection = 0;
	// Sampled parameter equal to the number of sites which have been removed from the data during ascertainment
	IntegerParameter ascSiteCount;
	//SnapSubstitutionModelGH m_substitutionmodel;
		
	
	@Override
	public void initAndValidate() {
		 
		speciesNetwork = speciesNetworkInput.get();
        Double sValues = m_pCoalescenceRate.get().getValue();
             
        String sCoalescenceRateValues = "";
        Double[] values = new Double[speciesNetwork.getBranchCount()];
    	 
        for (int i = 0; i < values.length; i++) {
                   	
            values[i] = new Double(sValues);                       
			sCoalescenceRateValues += values[i] + " ";
        }
                				 
        RealParameter coalescenceRate = m_pCoalescenceRate.get();
        coalescenceRate.initByName("value", sCoalescenceRateValues, "upper", 10.0, "lower", 0.0, "dimension", values.length);
		 
	}
		
	/**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
	 * @throws Exception 
     */
    @Override
    public double calculateLogP()  {
    	
	     
	    SanityChecks.checkNetworkSanity(speciesNetwork.getOrigin()); // species network should not be insane      
   
           
        // Handle snap data now !!!               
        m_pData = m_pDataInput.get();
       
        double u = m_pU.get().getValue();
        double v = m_pV.get().getValue();
                                    
        RealParameter  MyCoal = m_pCoalescenceRate.get();
        Double [] coalescenceRate = MyCoal.getValues();	 
       
        SnappNetLikelihoodCore m_core= new SnappNetLikelihoodCore(speciesNetwork,m_pData);
   
        
        try {
        	 
        		patternProbs=m_core.computeLogLikelihood(m_pData, speciesNetwork, u, v, coalescenceRate);       	        	
        		//compute likelihood for all sites summed up in patterns			 
			m_core=null;
 
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
         
        
        // calculate log prob      
     	logP = 0.0;
     	m_bUsenNonPolymorphic = m_usenNonPolymorphic.get();
   	
     	numPatterns = m_pData.getPatternCount();
     	for(int id = 0; id < numPatterns - (m_bUsenNonPolymorphic ? 0 : 2); id++) {
     		
     		double freq = m_pData.getPatternWeight(id);         		 
     		double siteL = patternProbs[id];
     		if (siteL==0.0) {
     			logP = -10e100;
     			break;
     		}
     		logP += (double)freq * Math.log(siteL);
     	}
      
     	
     	// correction for constant sites. If we are sampling the numbers of constant sites 
     	// (stored in ascSiteCount) then we include these probabilities. Otherwise we 
     	// assume that we want conditional likelihood, in which case we divide 
     	// by the probability that a site is not ascertained (or more correctly,
     	// subtract the log probability.
     	if (!m_bUsenNonPolymorphic) {
     				m_fP0 =  patternProbs[numPatterns - 2];
     				m_fP1 =  patternProbs[numPatterns - 1];
     				if (ascSiteCount != null) {   
     					ascLogP = (double)ascSiteCount.getValue(0) * Math.log(m_fP0) +
     							  (double)ascSiteCount.getValue(1) * Math.log(m_fP1);
     					logP += ascLogP;
     				} else {
     					logP -= (double) m_pData.getSiteCount() * Math.log(1.0 - m_fP0 - m_fP1);
     				}
     			}		
     	
 	
     	////////////////////////////////////
    	// calculate Likelihood Correction. 
		// When the assignment of individuals to populations/species is fixed, the allele counts in each population are sufficient 
		// statistics for the species tree parameters. However when testing species assignments this is no longer the case.
		// To address this we multiply the likelihood computed from allele counts by the probability of observing
		// the given sequences given those allele counts (and the species assignments).
		m_fLogLikelihoodCorrection = 0;
		if (useLogLikelihoodCorrection.get()) {
			// RRB: note that increasing the number of constant sites
			// does not change the m_fLogLikelihoodCorrection since the
			// contribution of constant sites is zero. This means,
			// m_fLogLikelihoodCorrection does not need to be recalculated
			// when ascSiteCount changes.
			// DJB: This is true, but only until we start looking at non-constant sites being ascertained.
	    	for (int i = 0; i < numPatterns; i++) {
	            int [] thisSite = m_pData.getPattern(i);  //count of red alleles for this site
	            int [] lineageCounts = m_pData.getPatternLineagCounts(i); //count of total lineages for this site
	            for (int j = 0; j < thisSite.length; j++) {
	            	m_fLogLikelihoodCorrection -= logBinom(thisSite[j], lineageCounts[j]) * m_pData.getPatternWeight(i);
	            }
	    	}
    	}		
		
	 	
		
		if (useLogLikelihoodCorrection.get()) {
			logP += m_fLogLikelihoodCorrection;
		}
		    		  	 
		
    	return logP;  	
    	
    }// calculateLogLikelihood

 
    private double logBinom(int k, int n) {
    	double f = 0;
    	for (int i = k + 1; i <= n; i++) {
    		f += Math.log(i) - Math.log(n - i + 1);
    	}
		return f;
	}

	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}


	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}




	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}
    
	
	


} // class SSSTreeLikelihood
