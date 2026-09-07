package Time::HiRes;

# Pure-Perl fallback for containers whose Perl lacks the compiled
# Time::HiRes (e.g. depot.galaxyproject.org's entrez-direct:22.4--he881be0_0,
# which ships perl-base without perl-modules -- confirmed 2026-09-07:
# `perl -MTime::HiRes -e1` fails with "Can't locate Time/HiRes.pm", which
# crashes every esearch/efetch call because edirect's own `nquire` helper
# uses Time::HiRes::time/usleep for NCBI rate-limit pacing. Implements only
# what edirect actually calls (confirmed via `grep -rho HiRes::`), plus the
# handful of common companions, via `date +%s.%N` and 4-arg select() --
# sub-second precision without an XS dependency. Not a general-purpose
# Time::HiRes replacement.
#
# Wired in via PERL5LIB pointing at this file's parent dir (see
# conf/provision_singularity.config's withLabel:'edirect' block); a real
# Time::HiRes already installed elsewhere on @INC always wins since PERL5LIB
# is searched in the order given and this shim is not injected there.

use strict;
use warnings;
use Exporter 'import';

our @EXPORT_OK = qw(time usleep sleep gettimeofday tv_interval);

sub time {
    my $t = `date +%s.%N`;
    chomp $t;
    return $t + 0 if $t;
    return CORE::time();
}

sub usleep {
    my ($usec) = @_;
    select(undef, undef, undef, $usec / 1_000_000);
    return $usec;
}

sub sleep {
    my ($sec) = @_;
    select(undef, undef, undef, $sec);
    return $sec;
}

sub gettimeofday {
    my $t = &time();
    my $sec = int($t);
    my $usec = int(($t - $sec) * 1_000_000);
    return wantarray ? ($sec, $usec) : $t;
}

sub tv_interval {
    my ($t0, $t1) = @_;
    $t1 ||= [gettimeofday()];
    return ($t1->[0] - $t0->[0]) + ($t1->[1] - $t0->[1]) / 1_000_000;
}

1;
