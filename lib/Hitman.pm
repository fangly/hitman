

=head1 OPTIONAL ARGUMENTS

=over

=item -l | --list

Simply list all available optional parameters.

=item -p <param> | --param <param>

Set an additional parameter. Re-use this argument for each parameter you want to
add, e.g. "-p SKIP_DENOISING=1 -p TRIM_LEN=250"

=for Euclid:
   repeatable
   param.type: string

=item -n <num_threads> | --num_threads <num_threads>

Specify the number of threads to use. Default: num_threads.default

=for Euclid:
   num_threads.type: +int
   num_threads.default: 12

=item -r | --report

Generate an HTML report.

=item -d <dir> | --dir <dir>

Specify the folder where to put the results. This folder will be created if
necessary. Default: dir.default

=for Euclid:
   dir.type: string
   dir.default: '.'

=back

=head1 VERSION

0.5

=head1 AUTHOR

Florent Angly <florent.angly@gmail.com>

=head1 COPYRIGHT

Copyright 2012-2014 Florent ANGLY <florent.angly@gmail.com>

Hitman is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Hitman is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Hitman. If not, see <http://www.gnu.org/licenses/>.

=cut


package Hitman;

use strict;
use warnings;
use File::Spec::Functions qw(catfile);
use Getopt::Euclid qw(:defer);
use IPC::System::Simple qw(runx);
use FindBin qw($RealBin $Script);


my $groovy_dir;
sub groovy_dir {
   # Get the location of the directory containing Bpipe Groovy files
   if (not defined $groovy_dir) {
      $groovy_dir = catfile($RealBin, '..', 'bpipe');
      if (not -d $groovy_dir) {
         die "Error: Directory '$groovy_dir' does not exist\n";
      }
   }
   return $groovy_dir;
}


my $groovy_script;
sub groovy_script {
   # Get the location of the Groovy file needed by Bpipe
   if (not defined $groovy_script) {
      $groovy_script = catfile(groovy_dir(), $Script);
      if (not -f $groovy_script) {
         die "Error: Groovy file '$groovy_script' does not exist\n";
      }
   }
   return $groovy_script;
}


sub parse_args {
   # Handle the special case where user wants to list available params
   my %args = map { $_ => undef } @ARGV;
   if (exists $args{'-l'} || exists $args{'--list'}) {
      list_params(groovy_script());
      exit;
   }
   # Let Getopt::Euclid handle the rest
   Getopt::Euclid->process_args(\@ARGV);
   return 1;
}


sub run_hitman {
   my ($inputs) = @_;
   # Prepare additional Bpipe parameters
   my @params;
   if (exists $ARGV{'--param'}) {
      for my $param (@{$ARGV{'--param'}}) {
         push @params, '-p', $param;
      }
   }

   # Set location of Bpipe shared modules
   $ENV{'BPIPE_LIB'} = groovy_dir();

   # Prepare Bpipe command, e.g.:
   #   bpipe run -p QUAL_TRUNC=13 -d out_dir pipeline.groovy infile1 infile2
   my @cmd = (
      'bpipe', 'run',
      @params,
      '-n', $ARGV{'--num_threads'},
      exists $ARGV{'--report'} ? '-r' : (),
      '-d', $ARGV{'--dir'},
      groovy_script(),
      @$inputs
   );

   # Run Bpipe
   runx( @cmd );

   ### TODO: DESTROY method to clean up:
   ###   bpipe stop
   ###   rm -rf .bpipe/ commandlog.txt
   ###   rm tmp_*

   return 1;
}


sub list_params {
   my ($groovy_script) = @_;
   # Read Groovy file and extract header section, between two stretch of
   # '///////' at the top of the file
   my $header;
   {
      local $/ = ('/'x80)."\n";
      open my $in, '<', $groovy_script or die "Error: Could not read file '$groovy_script': $!\n";
      <$in>;
      $header = <$in>;
      chomp $header;
      close $in;
   }
   # Parse header section and generate message.
   # Each parameter has a description line starting with '//' followed by
   # another line containing a key=value pair
   my $msg = '';
   while ( $header =~ m# ^ //\s*(.*?)\s* \n (\w+)\s*=\s*(.+?)\s* $ #mgx ) {
      my ($desc, $name, $default) = ($1, $2, $3);
      $msg .= "$name\n".
              "   $desc\n".
              "   Default: $default\n\n";
   }
   $msg ||= "No optional parameters available.\n";
   # Display message
   print $msg;
   return 1;
}


1;
