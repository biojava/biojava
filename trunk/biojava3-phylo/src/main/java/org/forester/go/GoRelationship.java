// $Id: GoRelationship.java,v 1.8 2008/12/06 01:26:51 cmzmasek Exp $
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

package org.forester.go;

public interface GoRelationship extends Comparable<GoRelationship> {

    public static final String PART_OF_STR              = "part_of";
    public static final String REGULATES_STR            = "regulates";
    public static final String NEGATIVELY_REGULATES_STR = "negatively_regulates";
    public static final String POSITIVELY_REGULATES_STR = "positively_regulates";

    public GoId getGoId();

    public Type getType();

    public static enum Type {
        PART_OF, REGULATES, NEGATIVELY_REGULATES, POSITIVELY_REGULATES;
    }
}
