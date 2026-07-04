;; Guix channels for reproducing the external-tool environment used by
;; `likegt geno` smoke tests.
;;
;; This pins the fetchable `guix-bioinformatics` state that still provides
;; `odgi` and `gafpack` on this system.  `gfainject` is checked separately by
;; scripts/check-geno-tools.sh because it is not available from the current
;; fetchable channel state.

(list
 (channel
  (name 'guix)
  (url "https://git.savannah.gnu.org/git/guix.git")
  (branch "master")
  (commit "44bbfc24e4bcc48d0e3343cd3d83452721af8c36"))
 (channel
  (name 'guix-bioinformatics)
  (url "https://gitlab.com/genenetwork/guix-bioinformatics.git")
  (branch "master")
  (commit "48af9393cf186230e08b0fa6f7f443bc818408d2")))
