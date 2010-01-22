// $Id: BranchWidth.java,v 1.9 2008/09/24 16:42:48 cmzmasek Exp $
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

package org.forester.phylogeny.data;

import java.io.IOException;
import java.io.Writer;

import org.forester.io.parsers.nhx.NHXtags;
import org.forester.io.parsers.phyloxml.PhyloXmlMapping;
import org.forester.util.ForesterUtil;

public class BranchWidth implements PhylogenyData {

    public final static double BRANCH_WIDTH_DEFAULT_VALUE = 1.0;
    private final double       _value;

    public BranchWidth() {
        _value = BRANCH_WIDTH_DEFAULT_VALUE;
    }

    public BranchWidth( final double value ) {
        _value = value;
    }

    @Override
    public StringBuffer asSimpleText() {
        return new StringBuffer( getValue() + "" );
    }

    @Override
    public StringBuffer asText() {
        return asSimpleText();
    }

    @Override
    public PhylogenyData copy() {
        return new BranchWidth( getValue() );
    }

    public double getValue() {
        return _value;
    }

    @Override
    public boolean isEqual( final PhylogenyData data ) {
        return getValue() == ( ( BranchWidth ) data ).getValue();
    }

    @Override
    public StringBuffer toNHX() {
        final StringBuffer sb = new StringBuffer();
        sb.append( NHXtags.PARENT_BRANCH_WIDTH );
        sb.append( getValue() );
        return sb;
    }

    @Override
    public void toPhyloXML( final Writer w, final int level, final String indentation ) throws IOException {
        w.write( ForesterUtil.LINE_SEPARATOR );
        w.write( indentation );
        PhylogenyDataUtil.appendElement( w, PhyloXmlMapping.WIDTH, getValue() + "" );
    }

    @Override
    public String toString() {
        return asText().toString();
    }
}
