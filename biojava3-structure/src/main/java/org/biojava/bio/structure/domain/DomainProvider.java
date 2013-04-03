package org.biojava.bio.structure.domain;

import java.util.SortedSet;

public interface DomainProvider {
	
	public SortedSet<String> getDomainNames(String name);
	
	public SortedSet<String> getRepresentativeDomains();
}
