// $Id: DomainSimilarityCalculator.java,v 1.10 2007/10/16 00:31:29 cmzmasek Exp
// $
//
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// All rights reserved
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA
//
// Contact: cmzmasek@yahoo.com
// WWW: www.phylosoft.org/forester

package org.forester.surfacing;

import java.util.List;
import java.util.SortedSet;

public interface DomainSimilarityCalculator {

    public SortedSet<DomainSimilarity> calculateSimilarities( final PairwiseDomainSimilarityCalculator pairwise_calculator,
                                                              final List<GenomeWideCombinableDomains> cdc_list,
                                                              final boolean ignore_domains_without_combinations_in_any_genome,
                                                              final boolean ignore_domains_specific_to_one_genome );;

    public static enum Detailedness {
        BASIC, LIST_COMBINING_DOMAIN_FOR_EACH_SPECIES, PUNCTILIOUS
    }

    public static enum GoAnnotationOutput {
        NONE, ALL
    }
}
