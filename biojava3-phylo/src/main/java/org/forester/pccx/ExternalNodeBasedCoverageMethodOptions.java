// $Id: ExternalNodeBasedCoverageMethodOptions.java,v 1.2 2008/03/09 00:11:19
// cmzmasek Exp $
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

package org.forester.pccx;

public class ExternalNodeBasedCoverageMethodOptions implements CoverageCalculationOptions {

    final private String _scoring_method;

    /**
     * This constructor sets the class name for the scoring method e.g.
     * "org.forester.tools.modeling.BranchCountingBasedScoringMethod"
     * 
     * @param scoring_method
     *            class name for the scoring method
     */
    public ExternalNodeBasedCoverageMethodOptions( final String scoring_method ) {
        _scoring_method = scoring_method;
    }

    public String asString() {
        final StringBuffer sb = new StringBuffer();
        sb.append( "scoring method: " );
        BranchCountingBasedScoringMethod scoring_method;
        try {
            scoring_method = ( BranchCountingBasedScoringMethod ) ( Class.forName( getScoringMethod() ) ).newInstance();
        }
        catch ( final Exception e ) {
            sb.append( "?" );
            return sb.toString();
        }
        sb.append( scoring_method.getDesciption() );
        return sb.toString();
    }

    public String getScoringMethod() {
        return _scoring_method;
    }
}
