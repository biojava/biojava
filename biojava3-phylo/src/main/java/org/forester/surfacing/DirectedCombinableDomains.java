// $Id: DirectedCombinableDomains.java,v 1.4 2008/08/20 01:15:36 cmzmasek Exp $
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

import java.util.ArrayList;
import java.util.List;

public class DirectedCombinableDomains extends BasicCombinableDomains {

    public DirectedCombinableDomains( final DomainId n_terminal_key_domain, final Species species ) {
        super( n_terminal_key_domain, species );
    }

    @Override
    public List<BinaryDomainCombination> toBinaryDomainCombinations() {
        final List<BinaryDomainCombination> binary_combinations = new ArrayList<BinaryDomainCombination>( getNumberOfCombinableDomains() );
        for( final DomainId domain : getCombiningDomains().keySet() ) {
            // Precondition (!): key domain is most upstream domain.
            //TODO ensure this is true.
            binary_combinations.add( new DirectedBinaryDomainCombination( getKeyDomain(), domain ) );
        }
        return binary_combinations;
    }
}
