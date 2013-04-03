// $Id: BranchData.java,v 1.16 2009/12/16 03:07:11 cmzmasek Exp $
// FORESTER -- software libraries and applications
// for evolutionary biology research and applications.
//
// Copyright (C) 2008-2009 Christian M. Zmasek
// Copyright (C) 2008-2009 Burnham Institute for Medical Research
// Copyright (C) 2000-2001 Washington University School of Medicine
// and Howard Hughes Medical Institute
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

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;

public class BranchData implements PhylogenyData {

    private BranchColor      _branch_color;
    private List<Confidence> _confidences;
    private BranchWidth      _branch_width;

    public BranchData() {
        // Doing nothing. 
    }

    public void addConfidence( final Confidence confidence ) {
        getConfidences().add( confidence );
    }

    @Override
    public StringBuffer asSimpleText() {
        throw new UnsupportedOperationException();
    }

    @Override
    public StringBuffer asText() {
        throw new UnsupportedOperationException();
    }

    @Override
    public PhylogenyData copy() {
        final BranchData new_bd = new BranchData();
        if ( isHasBranchColor() ) {
            new_bd.setBranchColor( ( BranchColor ) getBranchColor().copy() );
        }
        if ( isHasBranchWidth() ) {
            new_bd.setBranchWidth( ( BranchWidth ) getBranchWidth().copy() );
        }
        if ( isHasConfidences() ) {
            for( final Confidence confidence : getConfidences() ) {
                new_bd.addConfidence( ( Confidence ) confidence.copy() );
            }
        }
        return new_bd;
    }

    public BranchColor getBranchColor() {
        return _branch_color;
    }

    public BranchWidth getBranchWidth() {
        return _branch_width;
    }

    public Confidence getConfidence( final int index ) {
        return getConfidences().get( index );
    }

    public List<Confidence> getConfidences() {
        if ( _confidences == null ) {
            _confidences = new ArrayList<Confidence>();
        }
        return _confidences;
    }

    public int getNumberOfConfidences() {
        return getConfidences().size();
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        throw new UnsupportedOperationException();
    }

    public boolean isHasBranchColor() {
        return getBranchColor() != null;
    }

    public boolean isHasBranchWidth() {
        return getBranchWidth() != null;
    }

    public boolean isHasConfidences() {
        return getNumberOfConfidences() > 0;
    }

    public void setBranchColor( final BranchColor branch_color ) {
        _branch_color = branch_color;
    }

    public void setBranchWidth( final BranchWidth branch_width ) {
        _branch_width = branch_width;
    }

    @Override
    public StringBuffer toNHX() {
        final StringBuffer sb = new StringBuffer();
        if ( isHasConfidences() && ( getConfidence( 0 ).getValue() != Confidence.CONFIDENCE_DEFAULT_VALUE ) ) {
            sb.append( ":" );
            sb.append( getConfidence( 0 ).toNHX() );
        }
        if ( isHasBranchWidth() && ( getBranchWidth().getValue() != BranchWidth.BRANCH_WIDTH_DEFAULT_VALUE ) ) {
            sb.append( ":" );
            sb.append( getBranchWidth().toNHX() );
        }
        if ( isHasBranchColor() && ( getBranchColor().getValue() != null ) ) {
            sb.append( ":" );
            sb.append( getBranchColor().toNHX() );
        }
        return sb;
    }

    @Override
    public void toPhyloXML( final Writer writer, final int level, final String indentation ) throws IOException {
        if ( isHasConfidences() ) {
            for( final Confidence confidence : getConfidences() ) {
                confidence.toPhyloXML( writer, level, indentation );
            }
        }
        if ( isHasBranchWidth() ) {
            getBranchWidth().toPhyloXML( writer, level, indentation );
        }
        if ( isHasBranchColor() ) {
            getBranchColor().toPhyloXML( writer, level, indentation );
        }
    }
}
